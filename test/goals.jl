# Goals: the `with` / `without` keywords, the universe they are answered on, and
# the three-tier ladder of how much work an answer costs.

using Resolver: Bound, Dependency, Diagnosis, Fact, PkgData, Problem, Requirement,
    UserCompat, issatisfiable, pkg_info, prepare_pkg_info, resolve,
    resolve_prepared

GND() = Dict{Int,Vector{String}}()
GNC() = Dict{Int,Dict{String,Vector{Int}}}()

# A and B each need a different version of C, so they cannot coexist
const RIVALS = Dict(
    "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
    "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
    "C" => PkgData([2, 1], GND(), GNC()),
)

@testset "goals: the keyword shapes" begin
    # a package, a pair, and collections of either
    @test resolve(RIVALS, ["A"]; with = "C") == Dict("A" => 1, "C" => 1)
    @test resolve(RIVALS, ["A"]; with = "C" => 1) == Dict("A" => 1, "C" => 1)
    @test resolve(RIVALS, ["A"]; with = "C" => [1, 2]) == Dict("A" => 1, "C" => 1)
    @test resolve(RIVALS, ["A"]; with = ["A", "C"]) == Dict("A" => 1, "C" => 1)
    @test resolve(RIVALS, ["C"]; without = "A") == Dict("C" => 2)
    @test resolve(RIVALS, ["C"]; without = ["A", "B"]) == Dict("C" => 2)
    # a bare version reads as the singleton set containing it; a set is queried
    # with `in`, so an empty one admits nothing
    @test resolve(RIVALS, ["C"]; with = "C" => 2) == Dict("C" => 2)
    @test resolve(RIVALS, ["C"]; with = "C" => Int[], diagnose = false) ===
          nothing
    # both together
    @test resolve(RIVALS, ["C"]; with = "B", without = "A") ==
          Dict("B" => 1, "C" => 2)
    # no goal at all is the goal-⊤ path, untouched
    @test resolve(RIVALS, ["A"]; with = nothing, without = nothing) ==
          resolve(RIVALS, ["A"])

    # a `with` goal's package is in the answer, optimized like a requirement
    sol = resolve(RIVALS, String[]; with = "C")
    @test sol == Dict("C" => 2)
    # ... and a goal naming a package the universe does not have is
    # unsatisfiable, not silently dropped
    @test resolve(RIVALS, ["A"]; with = "Nonesuch", diagnose = false) === nothing
    @test resolve(RIVALS, ["A"]; without = "Nonesuch") == Dict("A" => 1, "C" => 1)
end

@testset "goals: a goal is never a fix" begin
    # the load-bearing sentence: `prob`'s constraints are negotiable, the goal
    # is the fixed ask. Adding B to `reqs` would make "drop requirement B" a
    # fix; asking for it with `with` must not.
    d = resolve(RIVALS, ["A"]; with = "B")
    @test d isa Diagnosis
    @test !isempty(d.fixes)
    @test all(a != Requirement("B") for f in d.fixes for a in f.actions)
    @test d.fixes == [d.fixes[1]]
    @test d.fixes[1].actions == Fact[Requirement("A")]
    # ... and the witness satisfies the goal, since the goal held over its own
    # resolve
    @test d.fixes[1].solution == Dict("B" => 1, "C" => 2)

    # the same problem with B as a requirement offers both directions
    d2 = resolve(RIVALS, ["A", "B"])
    @test Set(a for f in d2.fixes for a in f.actions) ==
          Set(Fact[Requirement("A"), Requirement("B")])

    # goal ≡ modified problem, on the verdict and on the solution -- the two
    # differ only in what a fix may propose
    for (reqs, w) in (("A", "C"), ("C", "A"), ("A", "B"))
        @test resolve(RIVALS, [reqs]; with = w, diagnose = false) ==
              resolve(RIVALS, sort([reqs, w]); diagnose = false)
    end
end

@testset "goals: the report names the goal" begin
    d = resolve(RIVALS, ["A"]; with = "B")
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("Conflict 1: the goal and A cannot be satisfied together.",
                   report)
    @test only(d.conflicts).goal
    @test only(d.conflicts).reqs == ["A"]
    @test occursin("A (all versions) works with C only at 1", report)
    @test occursin("B (all versions) works with C only at 2", report)

    # a conflict that has nothing to do with the goal says nothing about it
    d2 = resolve(RIVALS, ["A", "B"]; with = "C")
    @test !any(c.goal for c in d2.conflicts)
    @test occursin("Conflict 1: A, B cannot be satisfied together.",
                   sprint(show, MIME("text/plain"), d2))

    # a goal against the user's own constraints needs no requirement at all
    d3 = resolve(RIVALS, String[]; compat = Dict("C" => [1]), with = "B")
    @test only(d3.conflicts).goal
    @test isempty(only(d3.conflicts).reqs)
    report3 = sprint(show, MIME("text/plain"), d3)
    @test occursin("Conflict 1: the goal cannot be satisfied.", report3)
    @test occursin("your compat restricts C to 1", report3)
    @test occursin("relax your compat on C", report3)
end

@testset "goals: a `without` story is a chain of forced dependencies" begin
    # the Makie/FillArrays shape: nothing conflicts with anything, the banned
    # package is simply unavoidable, and the story is the forcing chain
    data = Dict(
        "M" => PkgData([1], Dict(1 => ["G"]), GNC()),
        "G" => PkgData([2, 1], Dict(2 => ["F"], 1 => ["F"]), GNC()),
        "F" => PkgData([1], GND(), GNC()),
        "P" => PkgData([1], GND(), GNC()),
    )
    d = resolve(data, ["M"]; without = "F")
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("Conflict 1: the goal and M cannot be satisfied together.",
                   report)
    @test occursin("M (all versions) requires G", report)
    @test occursin("G (all versions) requires F", report)
    @test Dependency("M", [1], "G", false) in only(d.conflicts).chain
    # G's two versions are interchangeable, so the universe carries one
    # representative of them -- the fact names what the instance has
    @test Dependency("G", [2], "F", false) in only(d.conflicts).chain
    @test d.fixes[1].actions == Fact[Requirement("M")]

    # an avoidable package is not a conflict: it is a different solution, and
    # the caller diffs the two to price it
    data2 = Dict(
        "M" => PkgData([2, 1], Dict(2 => ["F"]), GNC()),
        "F" => PkgData([1], GND(), GNC()),
    )
    @test resolve(data2, ["M"]) == Dict("M" => 2, "F" => 1)
    @test resolve(data2, ["M"]; without = "F") == Dict("M" => 1)

    # only *forced* edges are ever named: with a route around G the chain has
    # nothing honest to say, so it says nothing
    data3 = Dict(
        "M" => PkgData([2, 1], Dict(2 => ["G"]), GNC()),
        "G" => PkgData([1], Dict(1 => ["F"]), GNC()),
        "F" => PkgData([1], GND(), GNC()),
    )
    @test resolve(data3, ["M"]; without = "F") == Dict("M" => 1)
end

# The counterexample §4 of the design is built on: reachability keeps only the
# versions an *optimal* solution could use, which is exactly what a goal asks to
# step outside of. `A→B`, `B@3.1→X`, `B@2.3→∅`, no conflicts anywhere: the
# per-resolve filter keeps only B@3.1, so an instance prepared that way calls X
# unavoidable — while the raw problem avoids it by taking B@2.3.
const REACH = Dict(
    "A" => PkgData([1], Dict(1 => ["B"]), GNC()),
    "B" => PkgData([31, 23], Dict(31 => ["X"]), GNC()),
    "X" => PkgData([1], GND(), GNC()),
)

@testset "goals: the per-resolve instance answers the wrong question" begin
    prob = Problem(["A"])
    # the plain resolve takes the best B, and X comes along
    @test resolve(REACH, ["A"]) == Dict("A" => 1, "B" => 31, "X" => 1)

    # the per-resolve filter has dropped B@23 -- correctly, for its own
    # question: no optimal solution uses it
    prepared = prepare_pkg_info(pkg_info(REACH, prob), prob)
    @test prepared["B"].versions == [31]
    @test haskey(prepared, "X")
    # ... so asking *that* instance to avoid X is unsatisfiable, which is wrong
    @test resolve_prepared(prepared, Problem(["A", "X"]); diagnose = false) !==
          nothing
    sat = Resolver.SAT(prepared, prob, Resolver.Goal{String}(
        Pair{String,Any}[], ["X"]))
    try
        @test !Resolver.is_satisfiable(sat, ["A"])
    finally
        Resolver.finalize(sat)
    end

    # the goal route answers it: B@23 is in the T1 artifact, so avoiding X
    # costs a worse B and nothing else
    @test resolve(REACH, ["A"]; without = "X") == Dict("A" => 1, "B" => 23)
    @test issatisfiable(REACH, ["A"]; without = "X")
    # and the oracle agrees: a universe with X's only route deleted resolves
    # to exactly that
    free = Dict("A" => REACH["A"],
                "B" => PkgData([23], GND(), GNC()))
    @test resolve(free, ["A"]) == Dict("A" => 1, "B" => 23)

    # the same trap the other way round: a version goal below the reachable
    # prefix. Only B@31 survives preparation, so `B => 23` looks impossible
    @test resolve(REACH, ["A"]; with = "B" => 23) == Dict("A" => 1, "B" => 23)
    @test issatisfiable(REACH, ["A"]; with = "B" => 23)
end

# The brute-force reference for a goal: enumerate every valid solution, keep
# the ones the goal admits, and run the layered descent over what is left --
# `ref_resolve` with two filters and the goal's own packages among the descent's
# roots. A goal does not change what a solution *is*, only which ones count.
function ref_goal(
    data :: AbstractDict{P,<:PkgData{P,V}},
    reqs :: AbstractVector{P};
    with = Pair{P,Any}[],
    without = P[],
    by = identity,
) where {P,V}
    pkgs = sort!(collect(keys(data)))
    cands = Dict{P,V}[]
    TestResolver.each_potential_solution(data, pkgs) do s
        TestResolver.is_valid_solution(data, pkgs, s) || return
        push!(cands, Dict{P,V}(
            pkgs[i] => s[i] for i in eachindex(pkgs) if s[i] !== nothing))
    end
    filter!(c -> all(haskey(c, q) for q in reqs), cands)
    for (p, v) in with
        filter!(c -> haskey(c, p) && (v === nothing || c[p] == v), cands)
    end
    for p in without
        filter!(c -> !haskey(c, p), cands)
    end
    isempty(cands) && return nothing
    sol = Dict{P,V}()
    layer = sort!(unique(P[reqs; P[p for (p, _) in with]]); by)
    seen = Set{P}(layer)
    while true
        for p in layer
            rank = Dict{V,Int}(v => r for (r, v) in enumerate(data[p].versions))
            best = data[p].versions[minimum(rank[c[p]] for c in cands)]
            sol[p] = best
            filter!(c -> c[p] == best, cands)
        end
        next = P[]
        for p in layer
            haskey(data[p].depends, sol[p]) || continue
            for q in data[p].depends[sol[p]]
                q in seen || q in next || push!(next, q)
            end
        end
        isempty(next) && break
        layer = sort!(next; by)
        union!(seen, layer)
    end
    return sol
end

@testset "goals: brute force agrees on tiny instances" begin
    # every goal on every tiny random instance: `with`/`without` must agree
    # with resolving a universe that has the goal baked in by construction.
    Random.seed!(rand(RandomDevice(), UInt64))
    for (m, n) in ((2, 2), (3, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:150
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = [p for p = 1:m if rand(Bool)]
            isempty(reqs) && continue
            for p = 1:m
                # every goal shape against the enumerating oracle, and the
                # cheap predicate against both
                cases = Any[
                    ((without = p,), ref_goal(data, reqs; without = [p])),
                    ((with = p,),
                     ref_goal(data, reqs; with = Pair{Int,Any}[p => nothing])),
                ]
                for v in data[p].versions
                    push!(cases, ((with = p => v,),
                        ref_goal(data, reqs; with = Pair{Int,Any}[p => v])))
                end
                for (kws, ref) in cases
                    got = bare_resolve(data, reqs; kws...)
                    @test got == ref
                    @test (got !== nothing) == issatisfiable(data, reqs; kws...)
                end

                # ... and a `with` goal agrees with the modified problem, which
                # is the whole content of "a goal is not a constraint"
                @test bare_resolve(data, reqs; with = p) ==
                      bare_resolve(data, sort!(union(reqs, [p])))
            end
        end
    end
end

@testset "goals: the three-tier ladder agrees" begin
    # `issatisfiable`, `diagnose = false` and the full call answer the same
    # question with more and more work; nothing but the effort differs
    for (reqs, kws) in (
        (["A"], (;)),
        (["A"], (with = "B",)),
        (["A"], (with = "C" => 2,)),
        (["A"], (without = "C",)),
        (["C"], (without = "A",)),
        (["A", "B"], (;)),
    )
        yes = issatisfiable(RIVALS, reqs; kws...)
        probe = resolve(RIVALS, reqs; diagnose = false, kws...)
        full = resolve(RIVALS, reqs; kws...)
        @test yes == (probe !== nothing)
        @test yes == (full isa Dict)
        @test !yes || probe == full
        @test yes || full isa Diagnosis
    end

    # the probe really does skip the explanation
    @test resolve(RIVALS, ["A"]; with = "B", diagnose = false) === nothing

    # `issatisfiable` takes the same problem-or-requirements pair `resolve`
    # does, keywords included
    @test issatisfiable(RIVALS, Problem(["A"]))
    @test !issatisfiable(RIVALS, Problem(["A"]; compat = Dict("C" => [2])))
    @test !issatisfiable(RIVALS, ["A"]; compat = Dict("C" => [2]))
    @test !issatisfiable(RIVALS, ["A"]; pins = Dict("C" => 2))
end
