using Resolver
using Resolver: PkgData, diagnose, Diagnosis
using Test

# typed empty containers for P=String, V=Int, S=Vector{Int}
const NoDeps = Dict{Int,Vector{String}}
const NoComp = Dict{Int,Dict{String,Vector{Int}}}

# does resolve find a complete solution covering all requirements? resolve now
# returns a single solution covering all requirements, or nothing when they are
# not jointly satisfiable
resolve_complete(data, reqs) = resolve(data, reqs) !== nothing

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

@testset "solution pruning" begin
    NoUserCompat = Dict{String,Vector{Int}}
    NoHolds = Dict{String,Int}

    # SAT models may set unconstrained package variables true; reported
    # solutions must be pruned to the reachable closure of the satisfied
    # requirements. A/B fight over C (B via D); X has a dead bound on Z.
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["D"]), NoComp()),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
        "D" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "X" => PkgData([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => Int[]))),
        "Z" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    reqs = ["A", "B", "X"]
    d = diagnose(data, reqs)
    # first fix: greedy keeps A, drops B and X; solution is A's closure only —
    # no spurious D or Z
    fix = d.fixes[1]
    @test Set(a.pkg for a in fix.actions) == Set(["B", "X"])
    @test Set(keys(fix.solution)) == Set(["A", "C"])
    # every fix and upstream solution is closure-tight (checked in verifiers)
    for fix in d.fixes
        @test verify_fix(data, reqs, NoUserCompat(), NoHolds(), fix)
    end
    for cf in d.conflicts
        @test all(verify_upstream(data, cf.reqs, u) for u in cf.upstream)
    end
end

@testset "diagnose orchestration" begin
    # satisfiable input: nothing to diagnose
    sat_data = Dict("A" => PkgData([1], NoDeps(), NoComp()))
    @test_throws ArgumentError diagnose(sat_data, ["A"])

    # required package missing from data: throw, don't silently drop
    @test_throws ArgumentError("Required package B is not available in the package data") diagnose(sat_data, ["B"])

    # conflict caused solely by a stale third-party bound → upstream fix present
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => Int[]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    d = diagnose(data, ["A"])
    @test length(d.conflicts) == 1
    ups = d.conflicts[1].upstream
    @test !isempty(ups)
    @test all(verify_upstream(data, d.conflicts[1].reqs, u) for u in ups)

    # per-cluster probes: with a second, independent conflict present, the
    # bound-only conflict still gets its upstream suggestion
    data2 = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => Int[]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
        "X" => PkgData([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => Int[]))),
        "Z" => PkgData([1], NoDeps(), NoComp()),
    )
    d2 = diagnose(data2, ["A", "X"])
    @test length(d2.conflicts) == 2
    for cf in d2.conflicts
        @test !isempty(cf.upstream)
        @test all(verify_upstream(data2, cf.reqs, u) for u in cf.upstream)
    end

    # three-package chain: A→C∈{1}, B→D∈{3}, D@3→C∈{2,3}; C must be 1 and ≥2
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["D"]), Dict(1 => Dict("D" => [3]))),
        "C" => PkgData([1, 2, 3], NoDeps(), NoComp()),
        "D" => PkgData([1, 2, 3], Dict(3 => ["C"]), Dict(3 => Dict("C" => [2, 3]))),
    )
    d = diagnose(data, ["A", "B"])
    @test length(d.conflicts) == 1
    @test d.conflicts[1].reqs == ["A", "B"]
    chain_bounds = Resolver.Bound[f for f in d.conflicts[1].chain if f isa Resolver.Bound]
    @test chain_bounds == [
        Resolver.Bound("A", [1], "C", [1]),
        Resolver.Bound("B", [1], "D", [3]),
        Resolver.Bound("D", [3], "C", [2, 3]),
    ]
    @test length(d.fixes) >= 2
    for fix in d.fixes
        @test verify_fix(data, ["A", "B"],
            Dict{String,Vector{Int}}(), Dict{String,Int}(), fix)
    end
end

@testset "diagnose property (tiny)" begin
    for m = 1:2, n = 1:2
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            all(deps[i].bits >= deps[i+1].bits for i = 1:m-1) || continue
            for comp_bits = 0:2^c-1
                comp = make_comp(comp_bits)
                fill_data!(m, n, deps, comp, data)
                for reqs_bits = 1:2^m-1
                    reqs = collect(make_reqs(reqs_bits))
                    resolve_complete(data, reqs) && continue
                    diag = diagnose(data, reqs)
                    @test !isempty(diag.conflicts)
                    for cf in diag.conflicts
                        # cluster reqs are a requirement-level MUS
                        @test !resolve_complete(data, cf.reqs)
                        for p in cf.reqs
                            @test resolve_complete(data, filter(!=(p), cf.reqs))
                        end
                    end
                    for fix in diag.fixes
                        @test verify_fix(data, reqs,
                            Dict{Int,Vector{Int}}(), Dict{Int,Int}(), fix)
                    end
                end
            end
        end
    end
end

@testset "rendering" begin
    # format_versions run compression (V = Int)
    @test Resolver.format_versions([1, 2, 3, 4, 5], [1, 2, 3, 5]) == "1–3, 5"
    @test Resolver.format_versions([1, 2, 3], [1, 2, 3]) == "all versions"
    @test Resolver.format_versions([1, 2, 3], [2]) == "2"
    @test Resolver.format_versions([1, 2, 3], Int[]) == "no versions"
    @test Resolver.format_versions([10, 20, 30, 40], [10, 20, 40]) == "10–20, 40"
    # ranges read low–high (and runs ascend) even when `whole` is descending,
    # as the registry provider produces
    @test Resolver.format_versions([3, 2, 1], [1, 2]) == "1–2"
    @test Resolver.format_versions([5, 4, 3, 2, 1], [1, 2, 3, 5]) == "1–3, 5"

    # two-conflict example with a user compat, a bound, fixes and (per-cluster)
    # upstream suggestions for both conflicts
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1, 2]))),
        "C" => PkgData([1, 2, 3], NoDeps(), NoComp()),
        "X" => PkgData([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => Int[]))),
        "Z" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    d = diagnose(data, ["A", "X"]; compat = Dict("C" => [3]))
    @test length(d.conflicts) == 2
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("Conflict 1: A cannot be satisfied.", report)
    @test occursin("Conflict 2: X cannot be satisfied.", report)
    @test occursin("you require A", report)
    @test occursin("A (all versions) requires C at 1–2", report)
    @test occursin("your compat restricts C to 3", report)
    @test occursin("you require X", report)
    @test occursin("X (all versions) requires Z at no versions", report)
    @test occursin("Verified fixes:", report)
    @test occursin("relax your compat on C", report)
    @test occursin("drop requirement", report)
    @test occursin("→ allows: ", report)
    # the fix dropping every requirement allows nothing — rendered explicitly
    @test occursin("→ allows: nothing", report)
    @test !occursin("resolves:", report)
    @test !occursin("would give:", report)
    @test occursin("Upstream fixes:", report)
    @test occursin("a release of A relaxing its compat on C", report)
    @test occursin("a release of X relaxing its compat on Z", report)
    # compact summary (plural)
    @test sprint(show, d) == "Diagnosis(2 conflicts, $(length(d.fixes)) fixes)"

    # `allows` lines name only conflict-relevant packages (requirements +
    # contested dependencies), not the whole resolved closure; the complete
    # solution is still retained on the Fix struct. Here A and B need disjoint
    # versions of C, while Filler is an incidental dependency of A.
    data_trim = Dict(
        "A" => PkgData([1], Dict(1 => ["C", "Filler"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
        "Filler" => PkgData([1], NoDeps(), NoComp()),
    )
    dt = diagnose(data_trim, ["A", "B"])
    report_t = sprint(show, MIME("text/plain"), dt)
    @test occursin("→ allows: A 1, C 1", report_t)   # drop B
    @test occursin("→ allows: B 1, C 2", report_t)   # drop A
    @test !occursin("Filler", report_t)              # incidental dep omitted
    @test any(haskey(f.solution, "Filler") for f in dt.fixes)  # but kept on struct

    # single-conflict example with an upstream suggestion; also the singular
    # compact summary: one conflict, one fix (bound-only conflict)
    data1 = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), Dict(1 => Dict("B" => Int[]))),
        "B" => PkgData([1], NoDeps(), NoComp()),
    )
    d1 = diagnose(data1, ["A"])
    @test sprint(show, d1) == "Diagnosis(1 conflict, 1 fix)"
    report1 = sprint(show, MIME("text/plain"), d1)
    @test occursin("Conflict 1: A cannot be satisfied.", report1)
    @test occursin("A (all versions) requires B at no versions", report1)
    @test occursin("Verified fixes:", report1)
    @test occursin("drop requirement A", report1)
    @test occursin("Upstream fixes:", report1)
    @test occursin("a release of A relaxing its compat on B", report1)
end

end # @testset "diagnose"
