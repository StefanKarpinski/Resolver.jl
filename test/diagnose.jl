# UNSAT diagnostics.

using Resolver: Bound, Conflict, Diagnosis, Fact, Fix, Pin, PkgData, PkgInfo,
    Problem, Requirement, SAT, Uninstallable, UpstreamFix,
    UserCompat, dep_closure, filter_dominated, finalize, format_versions,
    Admission, MAX_UPSTREAM_SHOWN, interacting_pairs, order_chain, order_facts,
    pkg_info, prepare_pkg_info, render_action, restrict_info, superseded

@testset "diagnose: fact vocabulary" begin
    # facts compare and hash by value, across kinds
    @test Requirement("A") == Requirement("A")
    @test Requirement("A") ≠ Requirement("B")
    @test UserCompat("A", [1, 2]) == UserCompat("A", [1, 2])
    @test UserCompat("A", [1, 2]) ≠ UserCompat("A", [1])
    @test Pin("A", 1) == Pin("A", 1)
    @test Bound("A", [1], "C", [2]) == Bound("A", [1], "C", [2])
    @test Bound("A", [1], "C", [2]) ≠ Bound("A", [1], "C", [3])
    @test length(Set(Fact[Requirement("A"), Requirement("A")])) == 1
    @test length(Set(Fact[Requirement("A"), Uninstallable("A")])) == 2
    @test length(Set(Fact[Pin("A", 1), Pin("A", 2)])) == 2

    # deterministic ordering: actionability kind first, then package
    fs = Fact[Bound("B", [1], "C", [1]), Requirement("Z"), Pin("Y", 1),
              UserCompat("A", Int[]), Requirement("A"), Uninstallable("M")]
    @test order_facts(fs) == Fact[
        Requirement("A"), Requirement("Z"), Uninstallable("M"),
        Pin("Y", 1), UserCompat("A", Int[]), Bound("B", [1], "C", [1])]
end

@testset "diagnose: format_versions" begin
    @test format_versions([1, 2, 3, 4, 5], [1, 2, 3, 5]) == "1–3, 5"
    @test format_versions([1, 2, 3], [1, 2, 3]) == "all versions"
    @test format_versions([1, 2, 3], [2]) == "2"
    @test format_versions([1, 2, 3], Int[]) == "no versions"
    @test format_versions([10, 20, 30, 40], [10, 20, 40]) == "10–20, 40"
    # ranges read low–high (and runs ascend) even when `whole` is descending,
    # as the registry provider produces
    @test format_versions([3, 2, 1], [1, 2]) == "1–2"
    @test format_versions([5, 4, 3, 2, 1], [1, 2, 3, 5]) == "1–3, 5"
    # versions absent from `whole` are skipped, not an error: a fact may name a
    # version the filtered instance no longer carries
    @test format_versions([1, 2, 3], [2, 9]) == "2"
    @test format_versions([1, 2, 3], [9]) == "no versions"
    @test format_versions(Int[], [1]) == "no versions"
end

@testset "diagnose: filter_dominated" begin
    mkfix(actions...) = Fix{String,Int}(Fact[actions...], Dict{String,Int}())
    fA  = mkfix(Requirement("A"))
    fAB = mkfix(Requirement("A"), Requirement("B"))
    fBC = mkfix(Requirement("B"), UserCompat("C", Int[]))
    # a fix relaxing a proper subset of another's restrictions dominates it;
    # the dominated fix is dropped and the order of the rest preserved
    @test filter_dominated([fAB, fA, fBC]) == [fA, fBC]
    @test filter_dominated([fA, fBC]) == [fA, fBC]
    @test filter_dominated(Fix{String,Int}[]) == Fix{String,Int}[]
    # equal action sets don't dominate each other (⊊ is strict)
    @test length(filter_dominated([fA, mkfix(Requirement("A"))])) == 2
end

@testset "diagnose: order_chain" begin
    # story order: requirements first, then a BFS along the bound edges, each
    # package emitting its own facts as it is introduced
    facts = Fact[
        Bound("D", [3], "C", [2, 3]),
        Bound("B", [1], "D", [3]),
        Bound("A", [1], "C", [1]),
        Requirement("B"),
        Requirement("A"),
    ]
    chain = order_chain(facts, String, Int)
    @test chain == Fact[
        Requirement("A"), Requirement("B"),
        Bound("A", [1], "C", [1]),
        Bound("B", [1], "D", [3]),
        Bound("D", [3], "C", [2, 3]),
    ]

    # a package's own facts come out with it: pin then compat, before its bounds
    facts = Fact[
        UserCompat("C", [3]),
        Bound("A", [1], "C", [1, 2]),
        Requirement("A"),
        Pin("A", 1),
        Uninstallable("A"),
    ]
    chain = order_chain(facts, String, Int)
    @test chain == Fact[
        Requirement("A"), Uninstallable("A"), Pin("A", 1),
        Bound("A", [1], "C", [1, 2]), UserCompat("C", [3]),
    ]

    # facts the BFS can't reach are still emitted, in kind order
    facts = Fact[UserCompat("Z", Int[]), Requirement("A")]
    @test order_chain(facts, String, Int) ==
        Fact[Requirement("A"), UserCompat("Z", Int[])]

    # deterministic regardless of input order
    for _ = 1:8
        @test order_chain(shuffle(facts), String, Int) ==
            Fact[Requirement("A"), UserCompat("Z", Int[])]
    end
end

@testset "diagnose: rendering" begin
    versions = Dict("A" => [1], "C" => [1, 2, 3], "P" => [1, 2])
    d = Diagnosis{String,Int}(
        ["A"],
        [Conflict{String,Int}(["A"], Fact[
            Requirement("A"),
            Bound("A", [1], "C", [1, 2]),
            UserCompat("C", [3]),
            Pin("P", 2),
            Uninstallable("Q"),
        ], UpstreamFix{String,Int}[
            UpstreamFix(Bound("A", [1], "C", [1, 2]), Dict("A" => 1, "C" => 1)),
        ])],
        [Fix{String,Int}(Fact[UserCompat("C", [3])], Dict("A" => 1, "C" => 1)),
         Fix{String,Int}(Fact[Requirement("A")], Dict{String,Int}())],
        versions,
    )
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("Conflict 1: A cannot be satisfied.", report)
    @test occursin("you require A", report)
    @test occursin("A (all versions) works with C only at 1–2", report)
    @test occursin("your compat restricts C to 3", report)
    @test occursin("P is pinned at 2", report)
    @test occursin("no version of Q is installable", report)
    @test occursin("Verified fixes:", report)
    @test occursin("1. relax your compat on C", report)
    @test occursin("→ allows: A 1, C 1", report)
    @test occursin("2. drop requirement A", report)
    @test occursin("→ allows: nothing", report)
    @test occursin("Upstream fixes:", report)
    @test occursin("a release of A or C relaxing their compat on each other",
                   report)
    # requirements lead the "allows" list, then the other relevant packages
    @test !occursin("C 1, A 1", report)

    # compact summary: singular and plural
    @test sprint(show, d) == "Diagnosis(1 conflict, 2 fixes)"
    d1 = Diagnosis{String,Int}(["A"],
        [d.conflicts[1], d.conflicts[1]], [d.fixes[1]], versions)
    @test sprint(show, d1) == "Diagnosis(2 conflicts, 1 fix)"

    # a fix's solution is rendered only for conflict-relevant packages; the
    # complete assignment stays on the struct
    d2 = Diagnosis{String,Int}(["A"],
        [Conflict{String,Int}(["A"], Fact[Requirement("A")],
                              UpstreamFix{String,Int}[])],
        [Fix{String,Int}(Fact[Requirement("A")],
                         Dict("A" => 1, "Filler" => 1))],
        Dict("A" => [1], "Filler" => [1]))
    report2 = sprint(show, MIME("text/plain"), d2)
    @test occursin("→ allows: A 1", report2)
    @test !occursin("Filler", report2)
    @test d2.fixes[1].solution["Filler"] == 1

    # action rendering
    @test render_action(Requirement("A")) == "drop requirement A"
    @test render_action(UserCompat("A", Int[])) == "relax your compat on A"
    @test render_action(Pin("A", 1)) == "unpin A"
end

## diagnosis on the failed production instance

const NoDeps = Dict{Int,Vector{String}}
const NoComp = Dict{Int,Dict{String,Vector{Int}}}

# the diagnosis of an unsatisfiable resolve (and a check that it *is* one).
# `excludes` only travels in a `Problem` -- the convenience keywords carry
# compat and pins alone -- so this builds one
function diagnosis(data, reqs; excludes = [], kw...)
    d = resolve(data, Problem(reqs; excludes, kw...))
    @test d isa Diagnosis
    return d
end

# independently verify a fix against the *raw* problem: apply exactly the
# actions it names -- no more -- and check that the corrected problem really
# does resolve, to exactly the solution the fix recorded, covering every
# remaining requirement.
#
# "Exactly the actions it names" is the point. A fix relaxes a whole package's
# exclusion column, but only names the sources that put a cell in it; this
# checks that the ones it leaves unnamed really were doing nothing, so that
# following the report to the letter produces the promised versions.
function verify_fix(data, reqs, compat, pins, fix; excludes = [])
    reqs2 = collect(reqs)
    compat2 = Dict(compat)
    pins2 = Dict(pins)
    # (kind, package) pairs the fix says to relax
    offkind = Set{Tuple{Symbol,Any}}()
    for a in fix.actions
        a isa Requirement  && (reqs2 = filter(≠(a.pkg), reqs2))
        a isa UserCompat   && delete!(compat2, a.pkg)
        a isa Pin          && delete!(pins2, a.pkg)
        a isa Admission    && push!(offkind, (a.kind, a.pkg))
    end
    excludes2 = [Symbol(kind) => ((p, v) -> (Symbol(kind), p) ∉ offkind &&
                                            forbids(p, v)::Bool)
                 for (kind, forbids) in excludes]
    prob = Problem(reqs2; compat = compat2, pins = pins2, excludes = excludes2)
    sol = resolve(data, prob; diagnose = false)
    sol === nothing && return false
    sol == fix.solution && all(haskey(sol, p) for p in reqs2)
end

@testset "diagnose: requirement-level conflicts" begin
    # two requirements fighting over a shared dependency: one cluster, both
    # requirements in it, and dropping either one is a fix
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    d = diagnosis(data, ["A", "B"])
    @test length(d.conflicts) == 1
    @test d.conflicts[1].reqs == ["A", "B"]
    @test Set(Fact[Requirement("A"), Requirement("B")]) ⊆ Set(d.conflicts[1].chain)
    @test length(d.fixes) == 2
    @test Set(Set(a.pkg for a in f.actions) for f in d.fixes) ==
        Set([Set(["A"]), Set(["B"])])
    for fix in d.fixes
        @test verify_fix(data, ["A", "B"], Dict{String,Vector{Int}}(),
                         Dict{String,Int}(), fix)
    end

    # two independent conflicts: separate clusters, and every fix repairs both
    data2 = Dict(data...,
        "X" => PkgData([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => Int[]))),
        "Z" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    reqs = ["A", "B", "X"]
    d2 = diagnosis(data2, reqs)
    @test length(d2.conflicts) == 2
    @test [c.reqs for c in d2.conflicts] == [["A", "B"], ["X"]]
    @test !isempty(d2.fixes)
    for fix in d2.fixes
        @test any(a isa Requirement && a.pkg == "X" for a in fix.actions)
        @test verify_fix(data2, reqs, Dict{String,Vector{Int}}(),
                         Dict{String,Int}(), fix)
    end

    # cluster requirements are a requirement-level MUS: unsatisfiable, and
    # minimally so
    for c in d2.conflicts
        @test bare_resolve(data2, c.reqs) === nothing
        for p in c.reqs
            @test bare_resolve(data2, filter(≠(p), c.reqs)) !== nothing
        end
    end
end

@testset "diagnose: user constraints" begin
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), NoComp()),
        "B" => PkgData([1], NoDeps(), NoComp()),
    )
    compat = Dict("B" => Int[])
    d = diagnosis(data, ["A"]; compat)
    @test UserCompat("B", Int[]) in d.conflicts[1].chain
    # keep-order prefers keeping requirements, so relaxing the compat comes
    # before dropping the requirement
    @test length(d.fixes) == 2
    @test d.fixes[1].actions == Fact[UserCompat("B", Int[])]
    @test d.fixes[2].actions == Fact[Requirement("A")]
    @test d.fixes[1].solution == Dict("A" => 1, "B" => 1)
    @test isempty(d.fixes[2].solution)
    for fix in d.fixes
        @test verify_fix(data, ["A"], compat, Dict{String,Int}(), fix)
    end

    # a pin the user can lift
    d = diagnosis(data, ["A"]; pins = Dict("B" => 9))
    @test Pin("B", 9) in d.conflicts[1].chain
    @test d.fixes[1].actions == Fact[Pin("B", 9)]

    # a package carrying BOTH a compat bound and a pin is one relaxable group:
    # the two share an exclusion column, and relaxing half of one is exactly
    # the counterexample of Proposition D1'. The fix relaxes the whole column,
    # but names only the sources that put a cell in it -- and `verify_fix`
    # applies exactly what is named, so this is a test that the unnamed ones
    # really were inert.
    data3 = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), NoComp()),
        "B" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    compat3 = Dict("B" => [2])
    pins3 = Dict("B" => 1)
    d = diagnosis(data3, ["A"]; compat = compat3, pins = pins3)
    fix = d.fixes[1]
    @test Set(a.pkg for a in fix.actions) == Set(["B"])
    @test verify_fix(data3, ["A"], compat3, pins3, fix)

    # both sources named when both forbid something that is still there
    data4 = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), Dict(1 => Dict("B" => [2]))),
        "B" => PkgData([3, 2, 1], NoDeps(), NoComp()),
    )
    compat4 = Dict("B" => [1, 2])
    pins4 = Dict("B" => 1)
    d = diagnosis(data4, ["A"]; compat = compat4, pins = pins4)
    fix = d.fixes[1]
    @test any(a isa Pin for a in fix.actions)
    @test verify_fix(data4, ["A"], compat4, pins4, fix)

    # a constraint source that forbids nothing gets no selector, so it is not
    # named in the report either
    d = diagnosis(data3, ["A"];
                  compat = Dict("B" => [1, 2]), pins = Dict("B" => 9))
    @test !any(f isa UserCompat for f in d.conflicts[1].chain)
    @test any(f isa Pin for f in d.conflicts[1].chain)
end

@testset "diagnose: admission knobs" begin
    # An admission knob -- "no prereleases", "no yanked versions" -- is just
    # another constraint source, so it gets a selector, a fact and a fix like
    # any other. What it does not get is its own relaxation: it shares a
    # package's exclusion column with that package's compat bound and pin, and
    # Theorem D1 only licenses moving the column as a whole.
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), NoComp()),
        "B" => PkgData([2, 1], NoDeps(), NoComp()),
    )
    # every version of B is a "prerelease": nothing admissible is left. B's two
    # versions are interchangeable and forbidden alike, so the collapse keeps
    # one representative of them -- which is exactly why the fact names one
    # version rather than two
    nopre = ["prerelease" => (p, v) -> p == "B"] # a string kind is accepted
    d = diagnosis(data, ["A"]; excludes = nopre)
    chain = d.conflicts[1].chain
    @test Requirement("A") in chain
    @test Admission(:prerelease, "B", [2]) in chain
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("every version of B is a prerelease, which is not allowed",
                   report)
    @test occursin("allow prereleases of B", report)
    @test d.fixes[1].actions == Fact[Admission(:prerelease, "B", [2])]
    @test d.fixes[1].solution == Dict("A" => 1, "B" => 2)
    @test d.fixes[2].actions == Fact[Requirement("A")]

    # only some versions forbidden: the fact names them. A's own bound admits
    # only B@2, and B@2 is the prerelease -- so B@1 and B@3 survive the
    # collapse as one allowed class and B@2 as its own forbidden one
    data3 = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), Dict(1 => Dict("B" => [2]))),
        "B" => PkgData([3, 2, 1], NoDeps(), NoComp()),
    )
    some = Pair{Symbol,Any}[:prerelease => (p, v) -> p == "B" && v == 2]
    d = diagnosis(data3, ["A"]; excludes = some)
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("prereleases of B are not allowed (2)", report)
    @test d.fixes[1].actions == Fact[Admission(:prerelease, "B", [2])]
    @test d.fixes[1].solution == Dict("A" => 1, "B" => 2)
    @test verify_fix(data3, ["A"], Dict{String,Vector{Int}}(),
                     Dict{String,Int}(), d.fixes[1]; excludes = some)

    # a yanked newest version, with the older one ruled out by compat. Both
    # sources feed B's one exclusion column, but only the knob still forbids
    # anything the instance has -- the compat bound admits the survivor -- so
    # the fix names the knob alone, and `verify_fix` confirms that lifting just
    # the knob really does give the promised versions
    yanked = Pair{Symbol,Any}[:yanked => (p, v) -> p == "B" && v == 2]
    d = diagnosis(data, ["A"];
                  compat = Dict("B" => [2]), excludes = yanked)
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("every version of B is yanked", report)
    fix = d.fixes[1]
    @test fix.actions == Fact[Admission(:yanked, "B", [2])]
    @test occursin("allow the yanked version 2 of B", report)
    @test !occursin("relax your compat", report)
    @test fix.solution == Dict("A" => 1, "B" => 2)
    @test verify_fix(data, ["A"], Dict("B" => [2]), Dict{String,Int}(), fix;
                     excludes = yanked)

    # an unknown kind still renders sensibly rather than erroring
    odd = Pair{Symbol,Any}[:experimental => (p, v) -> p == "B"]
    d = diagnosis(data, ["A"]; excludes = odd)
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("every version of B is experimental, which is not allowed",
                   report)
    @test occursin("allow experimental versions of B", report)

    # a knob that forbids nothing anywhere leaves no trace
    inert = Pair{Symbol,Any}[:prerelease => (p, v) -> false]
    data2 = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    d = diagnosis(data2, ["A", "B"]; excludes = inert)
    @test !any(f isa Admission for f in d.conflicts[1].chain)
    @test !any(a isa Admission for f in d.fixes for a in f.actions)
end

@testset "diagnose: superseded prereleases" begin
    # the predicate: a prerelease is superseded by the release of the same base
    # version, and by nothing else
    vers = [v"1.2.3", v"1.2.3-alpha1", v"1.3.0-rc1", v"1.1.0"]
    @test  superseded(v"1.2.3-alpha1", vers)
    @test !superseded(v"1.3.0-rc1", vers)   # 1.3.0 is not out
    @test !superseded(v"1.2.3", vers)       # not a prerelease
    @test !superseded(v"1.2.3-alpha1", [v"1.2.3-alpha1", v"1.2.4"])
    # build metadata does not make a release a different version
    @test superseded(v"1.0.22-rc1", [v"1.0.22+0"])
    # a version type with no notion of prereleases is never superseded
    @test !superseded(2, [1, 2, 3])

    nodeps = Dict{VersionNumber,Vector{String}}()
    nocomp = Dict{VersionNumber,Dict{String,Vector{VersionNumber}}}()
    pre = Pair{Symbol,Any}[:prerelease => (p, v) -> !isempty(v.prerelease)]

    # not superseded: 2.0.0-rc1 is the only 2.0.0 there is, so admitting
    # prereleases is a real suggestion
    data = Dict(
        "A" => PkgData([v"1.0.0"], Dict(v"1.0.0" => ["B"]),
                       Dict(v"1.0.0" => Dict("B" => [v"2.0.0-rc1"]))),
        "B" => PkgData([v"2.0.0-rc1", v"1.0.0"], nodeps, nocomp),
    )
    d = diagnosis(data, ["A"]; excludes = pre)
    @test any(a isa Admission && a.kind === :prerelease
              for f in d.fixes for a in f.actions)
    @test occursin("allow prereleases of B",
                   sprint(show, MIME("text/plain"), d))

    # superseded: B@2.0.0 exists, so nobody should be told to install
    # 2.0.0-rc1 -- even though admitting it would technically resolve
    data2 = Dict(
        "A" => PkgData([v"1.0.0"], Dict(v"1.0.0" => ["B"]),
                       Dict(v"1.0.0" => Dict("B" => [v"2.0.0-rc1"]))),
        "B" => PkgData([v"2.0.0", v"2.0.0-rc1", v"1.0.0"], nodeps, nocomp),
    )
    d = diagnosis(data2, ["A"]; excludes = pre)
    report = sprint(show, MIME("text/plain"), d)
    @test !occursin("allow prereleases", report)
    @test !any(a isa Admission for f in d.fixes for a in f.actions)
end

@testset "diagnose: constraint labels" begin
    # "your compat" is wrong when the bound came from `Pkg.add(name, version)`
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), NoComp()),
        "B" => PkgData([2, 1], NoDeps(), NoComp()),
    )
    compat = Dict("B" => Int[])
    d = diagnosis(data, ["A"]; compat)
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("your compat restricts B to no versions", report)
    @test occursin("relax your compat on B", report)

    d = diagnosis(data, ["A"]; compat, labels = Dict("B" => :requested))
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("you requested B at no versions", report)
    @test occursin("relax the version you requested for B", report)
    @test !occursin("your compat", report)

    # an unknown label still renders, generically
    d = diagnosis(data, ["A"]; compat, labels = Dict("B" => Symbol("the lockfile")))
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("the lockfile restricts B to no versions", report)
    @test occursin("relax the the lockfile on B", report)
end

@testset "diagnose: uninstallable requirements" begin
    # a requirement the filter dropped outright -- no installable version --
    # cannot be reasoned about on the instance, so it gets its own conflict and
    # every fix has to drop it
    data = Dict(
        "A" => PkgData(Int[], NoDeps(), NoComp()),
        "B" => PkgData([1], NoDeps(), NoComp()),
    )
    d = diagnosis(data, ["A", "B"])
    @test length(d.conflicts) == 1
    @test d.conflicts[1].reqs == ["A"]
    @test d.conflicts[1].chain == Fact[Requirement("A"), Uninstallable("A")]
    @test length(d.fixes) == 1
    @test d.fixes[1].actions == Fact[Requirement("A")]
    @test d.fixes[1].solution == Dict("B" => 1)

    # transitively uninstallable, alongside a real conflict
    data2 = Dict(
        "A" => PkgData([1], Dict(1 => ["Gone"]), NoComp()),
        "Gone" => PkgData(Int[], NoDeps(), NoComp()),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => Int[]))),
        "C" => PkgData([1], NoDeps(), NoComp()),
    )
    d2 = diagnosis(data2, ["A", "B"])
    @test [c.reqs for c in d2.conflicts] == [["A"], ["B"]]
    @test !isempty(d2.fixes)
    for fix in d2.fixes
        @test Set(a.pkg for a in fix.actions) == Set(["A", "B"])
    end
end

@testset "diagnose: bound-level stories and upstream fixes" begin
    # A and B need disjoint versions of C: the story names both
    # incompatibilities, and relaxing either one alone would resolve it
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    d = diagnosis(data, ["A", "B"])
    bounds = Bound[f for f in d.conflicts[1].chain if f isa Bound]
    @test bounds == [Bound("A", [1], "C", [1]), Bound("B", [1], "C", [2])]
    ups = d.conflicts[1].upstream
    @test Set((u.bound.pkg, u.bound.dep) for u in ups) ==
        Set([("A", "C"), ("B", "C")])
    # every upstream suggestion carries a solution covering the cluster
    for u in ups
        @test all(haskey(u.solution, p) for p in d.conflicts[1].reqs)
    end

    # a conflict caused solely by a stale bound: one upstream suggestion
    data2 = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => Int[]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    d2 = diagnosis(data2, ["A"])
    @test length(d2.conflicts[1].upstream) == 1
    u = d2.conflicts[1].upstream[1]
    @test (u.bound.pkg, u.bound.dep) == ("A", "C")
    @test u.solution == Dict("A" => 1, "C" => 1)

    # a three-package chain: the story follows the dependency edges
    data3 = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["D"]), Dict(1 => Dict("D" => [3]))),
        "C" => PkgData([1, 2, 3], NoDeps(), NoComp()),
        "D" => PkgData([1, 2, 3], Dict(3 => ["C"]),
                       Dict(3 => Dict("C" => [2, 3]))),
    )
    d3 = diagnosis(data3, ["A", "B"])
    @test length(d3.conflicts) == 1
    bounds = Bound[f for f in d3.conflicts[1].chain if f isa Bound]
    # incompatibilities are undirected, so each is stated from its
    # alphabetically-first side: C@1 (which A needs) rules out the D@3 that B
    # needs, and D's own versions have collapsed to the two that differ
    @test bounds == [
        Bound("A", [1], "C", [1]),
        Bound("B", [1], "D", [3]),
        Bound("C", [1], "D", [1]),
    ]

    # bias: when the same conflict is tellable through a bound or through the
    # user's own compat, the requirement-level pass keeps the actionable one
    data4 = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => Int[]))),
        "C" => PkgData([1], NoDeps(), NoComp()),
    )
    d4 = diagnosis(data4, ["A"]; compat = Dict("C" => Int[]))
    @test UserCompat("C", Int[]) in d4.conflicts[1].chain
end

@testset "diagnose: incidental facts are demoted" begin
    # A minimal conflict has to close off the old versions of a package, which
    # means facts about versions the query already excludes. They belong in the
    # explanation but not at the same level as the line naming versions the
    # reader could actually get -- they are the story's tangents.
    #
    # B@3 is what the compat admits and it needs C@2; B@1 and B@2 are ruled out
    # by that same compat, so the facts closing them off are incidental.
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), NoComp()),
        "B" => PkgData([3, 2, 1], Dict(3 => ["C"], 2 => ["C"], 1 => ["C"]),
                       Dict(3 => Dict("C" => [2]),
                            2 => Dict("C" => [1]),
                            1 => Dict("C" => [1]))),
        "C" => PkgData([2, 1], NoDeps(), NoComp()),
    )
    compat = Dict("B" => [3], "C" => [1])
    d = diagnosis(data, ["A"]; compat)
    bounds = Bound[f for f in d.conflicts[1].chain if f isa Bound]
    @test length(bounds) == 2
    live = only(f for f in bounds if !f.incidental)
    @test live.versions == [3]          # the version the compat admits
    dead = only(f for f in bounds if f.incidental)
    # B@1 and B@2 say the same thing about C, so the collapse keeps one of them
    @test dead.versions == [2]          # a version the compat already excludes
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("B 3 works with C only at 2", report)
    @test occursin("(likewise for B 2, which your constraints already rule out)",
                   report)
    # the demoted fact gets no bullet of its own
    @test count("• B ", report) == 1

    # with *every* group incidental nothing is demoted: a verbose explanation
    # beats none at all
    d = diagnosis(data, ["A"]; compat = Dict("B" => Int[]))
    bounds = Bound[f for f in d.conflicts[1].chain if f isa Bound]
    @test all(!f.incidental for f in bounds)
end

@testset "diagnose: closure sub-instances" begin
    # the dependency closure unions over *all* versions' dependencies
    data = Dict(
        "A" => PkgData([1, 2], Dict(1 => ["B"], 2 => ["C"]), NoComp()),
        "B" => PkgData([1], NoDeps(), NoComp()),
        "C" => PkgData([1], NoDeps(), NoComp()),
        "Far" => PkgData([1], NoDeps(), NoComp()),
    )
    # (unfiltered, so that A@2 -- which reachability would drop, nothing
    # forcing A past its best version -- is still there to contribute C)
    info = pkg_info(data, ["A", "Far"]; filter = false)
    @test dep_closure(info, ["A"]) == Set(["A", "B", "C"])
    @test dep_closure(info, ["Far"]) == Set(["Far"])
    # restriction keeps the closure's structure intact
    sub = restrict_info(info, dep_closure(info, ["A"]))
    @test Set(keys(sub)) == Set(["A", "B", "C"])
    @test Resolver.check_info_structure(sub) === nothing
    @test bare_resolve(sub, ["A"]) == bare_resolve(info, ["A"])

    # conflicts with packages outside the closure are dropped, and by closure
    # exactness that changes no verdict
    data2 = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), NoComp()),
        "B" => PkgData([1], NoDeps(), NoComp()),
        "Out" => PkgData([1], NoDeps(), Dict(1 => Dict("B" => Int[]))),
    )
    info2 = pkg_info(data2, Problem(["A", "Out"]))
    sub2 = restrict_info(info2, dep_closure(info2, ["A"]))
    @test Set(keys(sub2)) == Set(["A", "B"])
    @test isempty(interacting_pairs(sub2))
    @test bare_resolve(sub2, ["A"]) == Dict("A" => 1, "B" => 1)
end

@testset "diagnose: the resolve API" begin
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([1, 2], NoDeps(), NoComp()),
    )
    prob = Problem(["A", "B"])
    info = pkg_info(data, prob) # a T1 artifact; `resolve` prepares it itself

    # every entry point returns the diagnosis itself, by default
    for d in (resolve(data, ["A", "B"]), resolve(data, prob),
              resolve(info, ["A", "B"]), resolve(info, prob))
        @test d isa Diagnosis{String,Int}
    end
    # ... and `nothing`, the bare verdict, when asked not to
    for d in (resolve(data, ["A", "B"]; diagnose = false),
              resolve(data, prob; diagnose = false),
              resolve(info, prob; diagnose = false))
        @test d === nothing
    end
    # show: a diagnosis at the REPL has to read as the answer it is
    d = resolve(data, prob)
    @test sprint(show, d) == "Diagnosis(1 conflict, 2 fixes)"
    report = sprint(show, MIME("text/plain"), d)
    @test startswith(report, "Unsatisfiable — 1 conflict, 2 fixes:\n")
    @test occursin("Conflict 1:", report)

    # a satisfiable resolve is unaffected
    @test resolve(data, ["A"]) == Dict("A" => 1, "C" => 1)

    # the frame discipline: diagnosing pops the selector frame and puts it
    # back, so a reused instance still resolves exactly as before
    cprob = Problem(["A", "B"]; compat = Dict("C" => [1]))
    # prepared exactly as `resolve` prepares it, so the instance under test is
    # the production one: collapsed to class representatives, then filtered
    cinfo = prepare_pkg_info(pkg_info(data, cprob), cprob)
    sat = SAT(cinfo, cprob)
    try
        before = resolve(sat, ["A"])
        @test before == Dict("A" => 1, "C" => 1)
        @test resolve(sat, ["A", "B"]) isa Diagnosis # diagnoses in place
        @test resolve(sat, ["A"]) == before          # ... and restores
        @test resolve(sat, ["B"]) isa Diagnosis      # C@2 still forbidden
        @test resolve(sat, ["A"]) == before
    finally
        finalize(sat)
    end
end

@testset "diagnose: property sweep on tiny data" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    for m = 2:3, n = 2:2
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        tested = 0
        for _ = 1:300
            tested < 40 || break
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            compat, pins = random_constraints(m, n)
            # an admission knob in the mix too, so the contract that a fix
            # names exactly the sources it needs is tested with all three
            excludes = rand(Bool) ? [] :
                let bad = Set{Int}(v for v = 1:n if rand(Bool))
                    Pair{Symbol,Any}[:test => (p, v) -> v ∈ bad]
                end
            reqs = collect(1:m)
            prob = Problem(reqs; compat, pins, excludes)
            resolve(data, prob; diagnose = false) === nothing || continue
            tested += 1
            dg = resolve(data, prob)
            @test !isempty(dg.conflicts)
            # every cluster is a requirement-level MUS of the raw problem
            for cf in dg.conflicts
                sub(rs) = Problem(rs; compat, pins, excludes)
                @test resolve(data, sub(cf.reqs); diagnose = false) === nothing
                for p in cf.reqs
                    @test resolve(data, sub(filter(≠(p), cf.reqs));
                                  diagnose = false) !== nothing
                end
            end
            # every fix really fixes it, to exactly the solution it reports,
            # applying exactly the actions it names and nothing more
            for fix in dg.fixes
                @test verify_fix(data, reqs, compat, pins, fix; excludes)
            end
            # no fix is dominated by another
            @test filter_dominated(dg.fixes) == dg.fixes
        end
        @test tested > 0
    end
end
