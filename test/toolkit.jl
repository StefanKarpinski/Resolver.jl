# The presentation toolkit: everything the default report decides, as something
# a caller can decide differently.

using Resolver: Bound, Change, Diagnosis, PkgData, Problem, Requirement,
    UserCompat, blame_phrase, changes, rank_upstream!, render_action,
    render_fact, report, resolve, superseded

# `sprint` does not forward keywords, and every knob here is one
rep(d; kws...) = sprint(io -> report(io, d; kws...))

TND() = Dict{Int,Vector{String}}()
TNC() = Dict{Int,Dict{String,Vector{Int}}}()

@testset "toolkit: show is report with the defaults" begin
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([2, 1], TND(), TNC()),
    )
    d = resolve(data, ["A", "B"])
    @test sprint(show, MIME("text/plain"), d) == rep(d)

    # max_upstream caps the list and says how many it dropped
    @test occursin("a release of", rep(d; max_upstream = 1))
    @test occursin("(1 more possible upstream fix not shown)",
                   rep(d; max_upstream = 1))
    @test !occursin("more possible upstream", rep(d))
    @test !occursin("not shown", rep(d; max_upstream = typemax(Int)))

    # trim_witnesses = false prints the whole assignment, not just the
    # contested packages
    full = rep(d; trim_witnesses = false)
    @test occursin("→ allows: A 1, C 1", rep(d))
    @test occursin("→ allows: A 1, C 1", full)  # nothing else here to add
end

@testset "toolkit: labels are a rendering choice too" begin
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), Dict(1 => Dict("B" => [2]))),
        "B" => PkgData([2, 1], TND(), TNC()),
    )
    prob = Problem(["A"]; compat = Dict("B" => [1]))
    d = resolve(data, prob)
    @test occursin("your compat restricts B to 1", rep(d))
    @test occursin("relax your compat on B", rep(d))
    # the same diagnosis, worded as `Pkg.add(name, version)` would want it
    said = rep(d; labels = Dict("B" => :requested))
    @test occursin("you requested B at 1", said)
    @test occursin("relax the version you requested for B", said)
    # ... and an unknown label reads neutrally rather than erroring
    @test occursin("lockfile restricts B to 1",
                   rep(d; labels = Dict("B" => :lockfile)))

    # `Problem`'s labels say the same thing at diagnosis time; the argument is
    # for callers whose wording is not a property of the problem
    d2 = resolve(data, Problem(["A"]; compat = Dict("B" => [1]),
                               labels = Dict("B" => :requested)))
    @test rep(d2) == said
end

@testset "toolkit: incidental facts, demoted or dropped" begin
    # a bound whose own side the user's compat already rules out: needed for
    # minimality, not what anyone is reading for
    data = Dict(
        "A" => PkgData([3, 2, 1], Dict(3 => ["B"], 2 => ["B"], 1 => ["B"]),
                       Dict(3 => Dict("B" => [2]), 2 => Dict("B" => Int[]),
                            1 => Dict("B" => Int[]))),
        "B" => PkgData([2, 1], TND(), TNC()),
    )
    prob = Problem(["A"]; compat = Dict("A" => [3], "B" => [1]))
    d = resolve(data, prob)
    @test d isa Diagnosis
    # whatever the demotion decides, the flag is on the fact for a client to
    # read, and every fact is still in the chain
    bounds = [f for c in d.conflicts for f in c.chain if f isa Bound]
    @test !isempty(bounds)
    plain = rep(d)
    verbose = rep(d; demote_incidental = false)
    @test length(split(verbose, "  • ")) ≥ length(split(plain, "  • "))
    @test !occursin("likewise for", verbose)
end

@testset "toolkit: the sentence layer" begin
    data = Dict(
        "A" => PkgData([1], Dict(1 => ["B"]), Dict(1 => Dict("B" => [2]))),
        "B" => PkgData([2, 1], TND(), TNC()),
    )
    d = resolve(data, Problem(["A"]; compat = Dict("B" => [1])))
    f = only(f for c in d.conflicts for f in c.chain if f isa UserCompat)
    @test sprint(render_fact, f, d) == "your compat restricts B to 1"
    @test render_action(f) == "relax your compat on B"
    @test blame_phrase(f) == "your compat on B"
    @test render_action(Requirement("A")) == "drop requirement A"
    @test blame_phrase(Requirement("A")) == "your requirement on A"
end

@testset "toolkit: superseded and upstream ranking" begin
    # the prerelease-supersession predicate, public so a client filtering its
    # own suggestions applies the same policy
    vers = [v"1.2.3", v"1.2.3-alpha1", v"1.3.0-beta1"]
    @test superseded(v"1.2.3-alpha1", vers)
    @test !superseded(v"1.3.0-beta1", vers)
    @test !superseded(v"1.2.3", vers)
    @test !superseded(1, [1, 2])   # a version type with no prerelease notion

    data = Dict(
        "A" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
        "B" => PkgData([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
        "C" => PkgData([2, 1], TND(), TNC()),
    )
    d = resolve(data, ["A", "B"])
    ups = [u for c in d.conflicts for u in c.upstream]
    @test length(ups) == 2
    # a client's own filter, over public structure
    mine = filter(u -> u.bound.pkg == "A" || u.bound.dep == "A", ups)
    @test !isempty(mine)
    # ranking off the diagnosis's own version lists, which is all a client
    # holding a report has
    ranked = rank_upstream!(copy(ups), d)
    @test Set(ranked) == Set(ups)
    @test rank_upstream!(copy(ranked), d) == ranked   # deterministic
end

@testset "toolkit: changes prices a feasible goal" begin
    # the `without = "TranscodingStreams"` shape: not a conflict, a different
    # solution, and the caller diffs the two to say what it costs
    data = Dict(
        "CSV" => PkgData([10, 5], Dict(10 => ["Codec"]), TNC()),
        "Codec" => PkgData([1], TND(), TNC()),
    )
    plain = resolve(data, ["CSV"])
    free = resolve(data, ["CSV"]; without = "Codec")
    @test plain == Dict("CSV" => 10, "Codec" => 1)
    @test free == Dict("CSV" => 5)
    cs = changes(plain, free)
    @test cs == [Change("CSV", 10, 5), Change("Codec", 1, nothing)]
    @test sprint(show, cs[1]) == "CSV 10 → 5"
    @test sprint(show, cs[2]) == "Codec 1 removed"
    @test sprint(show, Change("X", nothing, 2)) == "X 2 added"
    # a diff with itself is empty, and the diff is antisymmetric
    @test isempty(changes(plain, plain))
    @test changes(free, plain) == [Change("CSV", 5, 10), Change("Codec", nothing, 1)]
end
