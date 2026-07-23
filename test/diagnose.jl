using Resolver
using Resolver: PkgData, Diagnosis
using Test

# typed empty containers for P=String, V=Int, S=Vector{Int}
const NoDeps = Dict{Int,Vector{String}}
const NoComp = Dict{Int,Dict{String,Vector{Int}}}

# does resolve find a complete solution covering all requirements?
function resolve_complete(data, reqs)
    pkgs, vers = resolve(data, reqs)
    size(vers, 2) == 0 && return isempty(reqs)
    ridx = indexin(collect(reqs), pkgs)
    any(all(!isnothing(vers[i, k]) for i in ridx) for k in axes(vers, 2))
end

# build a diagnostic instance, run body, free it
function with_diag(body, data, reqs; kw...)
    diag = Resolver.DiagSAT(data, reqs; kw...)
    try body(diag)
    finally
        Resolver.finalize(diag)
    end
end

# is the whole diagnostic instance (all selectors on) satisfiable?
diag_sat_all(data, reqs; kw...) =
    with_diag(data, reqs; kw...) do diag
        Resolver.sat_diag(diag, Resolver.all_selectors(diag))
    end

# is a solution dict valid for the given data?
function valid_solution(data, sol)
    pkgs = collect(keys(sol))
    vers = Union{Int,Nothing}[sol[p] for p in pkgs]
    TestResolver.is_valid_solution(data, pkgs, vers)
end

# is a solution closure-tight? i.e. every package in it is dependency-reachable
# from the given roots along the chosen versions' depends edges
function closure_tight(data, roots, sol)
    seen = Set{keytype(sol)}(p for p in roots if haskey(sol, p))
    queue = collect(seen)
    while !isempty(queue)
        p = pop!(queue)
        haskey(data[p].depends, sol[p]) || continue
        for q in data[p].depends[sol[p]]
            haskey(sol, q) && q ∉ seen || continue
            push!(seen, q)
            push!(queue, q)
        end
    end
    Set(keys(sol)) == seen
end

# independently verify a fix: apply its actions (drop reqs / relax compat /
# unhold), check the modified problem is satisfiable, and check the fix's
# recorded solution is valid, respects residual user constraints and covers all
# remaining requirements.
function verify_fix(data, reqs, compat, holds, fix)
    reqs2 = collect(reqs)
    compat2 = Dict(compat)
    holds2 = Dict(holds)
    for a in fix.actions
        if a isa Resolver.Requirement
            reqs2 = filter(!=(a.pkg), reqs2)
        elseif a isa Resolver.UserCompat
            delete!(compat2, a.pkg)
        elseif a isa Resolver.Hold
            delete!(holds2, a.pkg)
        end
    end
    ok = diag_sat_all(data, reqs2; compat = compat2, holds = holds2)
    sol = fix.solution
    valid = valid_solution(data, sol)
    for (p, set) in compat2
        haskey(sol, p) && !(sol[p] in set) && (valid = false)
    end
    for (p, v) in holds2
        haskey(sol, p) && sol[p] != v && (valid = false)
    end
    covers = all(haskey(sol, p) for p in reqs2)
    tight = closure_tight(data, reqs2, sol)
    return ok && valid && covers && tight
end

# verify an upstream fix: relaxing the named bound (dropping pkg's compat on dep
# for the sharing versions) makes this conflict's cluster satisfiable, and the
# recorded solution is valid against that relaxed data and covers the cluster.
function verify_upstream(data, creqs, u)
    data2 = deepcopy(data)
    dp = data2[u.bound.pkg]
    for v in u.bound.versions
        haskey(dp.compat, v) && delete!(dp.compat[v], u.bound.dep)
    end
    ok = diag_sat_all(data2, creqs)
    valid = valid_solution(data2, u.solution)
    covers = all(haskey(u.solution, p) for p in creqs)
    tight = closure_tight(data2, creqs, u.solution)
    return ok && valid && covers && tight
end

@testset "diagnose" begin

@testset "builder equivalence" begin
    # a batch of no-compat/no-holds cases: diagnostic all-selectors-on verdict
    # must match whether resolve finds a complete solution.
    cases = []

    # 1. simple SAT: A needs C, C available
    push!(cases, (
        Dict(
            "A" => PkgData([1], Dict(1 => ["C"]), NoComp()),
            "C" => PkgData([1, 2], NoDeps(), NoComp()),
        ), ["A"]))

    # 2. SAT via satisfiable third-party bound: A@1 allows C∈{2}
    push!(cases, (
        Dict(
            "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
            "C" => PkgData([1, 2], NoDeps(), NoComp()),
        ), ["A"]))

    # 3. UNSAT: A and B impose conflicting bounds on shared dep C
    push!(cases, (
        Dict(
            "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
            "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
            "C" => PkgData([1, 2], NoDeps(), NoComp()),
        ), ["A", "B"]))

    # 4. UNSAT reverse direction: C requires A but forbids every A version
    push!(cases, (
        Dict(
            "A" => PkgData([1, 2], NoDeps(), NoComp()),
            "C" => PkgData([1], Dict(1 => ["A"]), Dict(1 => Dict("A" => Int[]))),
        ), ["C"]))

    # 5. SAT reverse direction: C requires A, C@1 allows only A∈{2}
    push!(cases, (
        Dict(
            "A" => PkgData([1, 2], NoDeps(), NoComp()),
            "C" => PkgData([1], Dict(1 => ["A"]), Dict(1 => Dict("A" => [2]))),
        ), ["C"]))

    for (data, reqs) in cases
        @test diag_sat_all(data, reqs) == resolve_complete(data, reqs)
    end

    # user compat can turn a satisfiable problem unsatisfiable
    data = Dict(
        "A" => PkgData([1, 2], Dict(1 => ["C"], 2 => ["C"]), NoComp()),
        "C" => PkgData([1], NoDeps(), NoComp()),
    )
    @test resolve_complete(data, ["A"])
    @test diag_sat_all(data, ["A"])
    # forbid every version of A via user compat
    @test !diag_sat_all(data, ["A"]; compat = Dict("A" => Int[]))

    # a manifest hold can too: hold C to a nonexistent version
    @test !diag_sat_all(data, ["A"]; holds = Dict("C" => 9))

    # a hold at a package's only version excludes nothing → no H group
    with_diag(data, ["A"]; holds = Dict("C" => 1)) do diag
        @test isempty(diag.Hs)
    end
end

@testset "B coalescing" begin
    # P has 4 versions: v1-2 share one bound on Q, v3-4 share another → 2 groups
    data = Dict(
        "P" => PkgData([1, 2, 3, 4], NoDeps(), Dict(
            1 => Dict("Q" => [1, 2]),
            2 => Dict("Q" => [1, 2]),
            3 => Dict("Q" => [2, 3]),
            4 => Dict("Q" => [2, 3]),
        )),
        "Q" => PkgData([1, 2, 3], NoDeps(), NoComp()),
    )
    with_diag(data, ["P"]) do diag
        n = count(f isa Resolver.Bound && f.pkg == "P" && f.dep == "Q"
                  for (_, f) in diag.Bs)
        @test n == 2
    end

    # versions with no compat entry, or a compat entry that excludes nothing,
    # join no group → only 1 group here
    data = Dict(
        "P" => PkgData([1, 2, 3], NoDeps(), Dict(
            1 => Dict("Q" => [1, 2]),      # excludes {3} → group
            3 => Dict("Q" => [1, 2, 3]),   # excludes nothing → no group
        )),                                # version 2 has no entry → no group
        "Q" => PkgData([1, 2, 3], NoDeps(), NoComp()),
    )
    with_diag(data, ["P"]) do diag
        n = count(f isa Resolver.Bound && f.pkg == "P" && f.dep == "Q"
                  for (_, f) in diag.Bs)
        @test n == 1
    end
end

@testset "group MUS" begin
    # direct conflict: A requires C∈{1}, B requires C∈{2}
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    with_diag(data, ["A", "B"]) do diag
        sels = Resolver.all_selectors(diag)
        mus = Resolver.group_mus(diag, sels, Resolver.full_shrink(diag))
        facts = Set(Resolver.facts_of(diag, mus))
        @test facts == Set(Resolver.Fact[
            Resolver.Requirement("A"),
            Resolver.Requirement("B"),
            Resolver.Bound("A", [1], "C", [1]),
            Resolver.Bound("B", [1], "C", [2]),
        ])
        # minimality: dropping any single member is satisfiable
        for s in mus
            @test Resolver.sat_diag(diag, filter(!=(s), mus))
        end
    end

    # bias: the same conflict is tellable through a bound OR a user compat entry
    # (both forbid the sole version of C); the MUS keeps the actionable compat.
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => Int[]))),
        "C" => PkgData([1], NoDeps(), NoComp()),
    )
    with_diag(data, ["A"]; compat = Dict("C" => Int[])) do diag
        sels = Resolver.all_selectors(diag)
        mus = Resolver.group_mus(diag, sels, Resolver.full_shrink(diag))
        facts = Set(Resolver.facts_of(diag, mus))
        @test facts == Set(Resolver.Fact[
            Resolver.Requirement("A"),
            Resolver.UserCompat("C", Int[]),
        ])
        @test !any(f isa Resolver.Bound for f in facts)
    end
end

@testset "requirement clustering" begin
    # two independent conflicts on disjoint packages:
    #   A vs B fight over shared dep C; X vs Y fight over shared dep Z
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
        "X" => PkgData([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => [1]))),
        "Y" => PkgData([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => [2]))),
        "Z" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    with_diag(data, ["A", "B", "X", "Y"]) do diag
        clusters = Resolver.cluster_reqs(diag)
        pkgsets = Set(Set(diag.facts[s].pkg for s in c) for c in clusters)
        @test pkgsets == Set([Set(["A", "B"]), Set(["X", "Y"])])
    end

    # singleton cluster: A is individually uninstallable because the user's
    # compat excludes every version of it
    data = Dict("A" => PkgData([1, 2], NoDeps(), NoComp()))
    with_diag(data, ["A"]; compat = Dict("A" => Int[])) do diag
        clusters = Resolver.cluster_reqs(diag)
        @test length(clusters) == 1
        @test Set(diag.facts[s].pkg for s in clusters[1]) == Set(["A"])
        # its group MUS is exactly {R(A), C(A)}
        cluster = clusters[1]
        sels_on = Int[cluster; Resolver.compat_sels(diag);
                      Resolver.hold_sels(diag); Resolver.bound_sels(diag)]
        mus = Resolver.group_mus(diag, sels_on, Resolver.full_shrink(diag, cluster))
        facts = Set(Resolver.facts_of(diag, mus))
        @test facts == Set(Resolver.Fact[
            Resolver.Requirement("A"),
            Resolver.UserCompat("A", Int[]),
        ])
    end
end

@testset "fix enumeration" begin
    NoUserCompat = Dict{String,Vector{Int}}
    NoHolds = Dict{String,Int}

    # preference: both "relax compat on B" and "drop requirement A" fix the
    # problem; keep-order keeps requirements and drops compat first.
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), NoComp()),
        "B" => PkgData([1], NoDeps(), NoComp()),
    )
    compat = Dict("B" => Int[])
    with_diag(data, ["A"]; compat) do diag
        fixes = Resolver.enumerate_fixes(diag)
        @test length(fixes) >= 2
        @test length(fixes[1].actions) == 1
        @test fixes[1].actions[1] isa Resolver.UserCompat
        @test fixes[1].actions[1].pkg == "B"
        @test any(a isa Resolver.Requirement && a.pkg == "A"
                  for a in fixes[2].actions)
        for fix in fixes
            @test verify_fix(data, ["A"], compat, NoHolds(), fix)
        end
    end

    # max_fixes respected
    with_diag(data, ["A"]; compat) do diag
        @test length(Resolver.enumerate_fixes(diag; max_fixes = 1)) == 1
    end

    # multi-cluster: two independent conflicts; every fix must repair BOTH
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
        "X" => PkgData([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => [1]))),
        "Y" => PkgData([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => [2]))),
        "Z" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    reqs = ["A", "B", "X", "Y"]
    with_diag(data, reqs) do diag
        fixes = Resolver.enumerate_fixes(diag)
        @test !isempty(fixes)
        for fix in fixes
            @test verify_fix(data, reqs, NoUserCompat(), NoHolds(), fix)
        end
    end

    # dominated-fix suppression: a fix removing a proper subset of another's
    # restrictions dominates it; the dominated fix is dropped, order preserved
    mkfix(actions...) = Resolver.Fix{String,Int}(
        Resolver.Fact[actions...], Dict{String,Int}())
    fA  = mkfix(Resolver.Requirement("A"))
    fAB = mkfix(Resolver.Requirement("A"), Resolver.Requirement("B"))
    fBC = mkfix(Resolver.Requirement("B"), Resolver.UserCompat("C", Int[]))
    @test Resolver.filter_dominated([fAB, fA, fBC]) == [fA, fBC]
    @test Resolver.filter_dominated([fA, fBC]) == [fA, fBC]

    # no fix returned by enumerate_fixes is dominated by another
    with_diag(data, reqs) do diag
        fixes = Resolver.enumerate_fixes(diag)
        @test Resolver.filter_dominated(fixes) == fixes
    end
end

end # @testset "diagnose"
