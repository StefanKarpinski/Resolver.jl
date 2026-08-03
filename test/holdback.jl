# Diagnostics for resolutions that succeeded and still surprised someone.

using Resolver: Admission, Bound, Fact, Holdback, Holdbacks, Pin, PkgData,
    Problem, Requirement, UserCompat, holdbacks, improved, pkg_info,
    prepare_pkg_info, resolve_prepared, summarize

const HND = Dict{Int,Vector{String}}
const HNC = Dict{Int,Dict{String,Vector{Int}}}

# the shape of the surprise this exists for: a "julia" whose every release
# bundles exactly one version of a stdlib, so a bound on the stdlib steers the
# whole toolchain
const STEERING = Dict(
    "App"   => PkgData([1], Dict(1 => ["julia"]), HNC()),
    "julia" => PkgData([12, 11, 10],
                       Dict(12 => ["Stats"], 11 => ["Stats"], 10 => ["Stats"]),
                       Dict(12 => Dict("Stats" => [12]),
                            11 => Dict("Stats" => [11]),
                            10 => Dict("Stats" => [10]))),
    "Stats" => PkgData([12, 11, 10], HND(), HNC()),
)

# resolve, then explain -- on the same prepared universe the solution came from
function held(data, prob, pkgs = nothing; kw...)
    info = prepare_pkg_info(pkg_info(data, prob), prob)
    sol = resolve_prepared(info, prob)
    @test sol !== nothing
    hs = pkgs === nothing ? holdbacks(info, prob, sol; kw...) :
                            holdbacks(info, prob, sol, pkgs; kw...)
    return sol, hs
end

@testset "holdback: a stale bound steers the toolchain" begin
    prob = Problem(["App", "julia"]; compat = Dict("Stats" => [10]))
    sol, hs = held(STEERING, prob, ["julia"])
    @test sol["julia"] == 10
    h = only(hs)
    @test h.pkg == "julia"
    @test h.resolved == 10
    @test h.best == 12
    # the story names the bound and the incompatibility it creates
    @test UserCompat("Stats", [10]) in h.chain
    @test Bound("julia", [12], "Stats", [12]) in h.chain
    # the fix is verified, and its witness is the optimal resolution under it
    @test h.fix !== nothing
    @test h.fix.actions == Fact[UserCompat("Stats", [10])]
    @test improved(h) == 12
    @test h.fix.solution == Dict("App" => 1, "julia" => 12, "Stats" => 12)
    relaxed = Problem(["App", "julia"])
    @test h.fix.solution ==
        resolve_prepared(prepare_pkg_info(pkg_info(STEERING, relaxed), relaxed),
                         relaxed)

    report = sprint(show, MIME("text/plain"), h)
    @test occursin("julia resolved to 10; 12 is available.", report)
    @test occursin("your compat restricts Stats to 10", report)
    @test occursin("julia 12 works with Stats only at 12", report)
    @test occursin("⇒ held back by your compat on Stats.", report)
    @test occursin("relax your compat on Stats → allows: julia 12, Stats 12",
                   report)
    @test summarize(h) ==
        "julia is held back to 10 by your compat on Stats; relaxing it allows julia 12"
    @test sprint(show, h) == "Holdback(julia 10 < 12, fixable)"
end

@testset "holdback: nothing to explain" begin
    # unconstrained: everything is at its best, so there is no holdback and no
    # solving done to discover it
    prob = Problem(["App", "julia"])
    sol, hs = held(STEERING, prob)
    @test sol["julia"] == 12
    @test isempty(hs)

    # a package that could have been at its best but for the priority order is
    # not "held back" by anything a report could name
    data = Dict(
        "A" => PkgData([2, 1], HND(), Dict(2 => Dict("B" => [1]))),
        "B" => PkgData([2, 1], HND(), HNC()),
    )
    prob = Problem(["A", "B"])
    sol, hs = held(data, prob)
    @test sol == Dict("A" => 2, "B" => 1)   # A wins, B settles
    @test isempty(hs)                       # ... and no constraint to blame
end

@testset "holdback: every kind of constraint" begin
    # a pin
    prob = Problem(["App", "julia"]; pins = Dict("Stats" => 11))
    _, hs = held(STEERING, prob, ["julia"])
    h = only(hs)
    @test h.resolved == 11
    @test Pin("Stats", 11) in h.chain
    @test h.fix.actions == Fact[Pin("Stats", 11)]
    @test improved(h) == 12
    @test occursin("⇒ held back by your pin on Stats.",
                   sprint(show, MIME("text/plain"), h))

    # a relaxable admission kind
    prob = Problem(["App", "julia"];
                   excludes = [:prerelease => (p, v) -> p == "Stats" && v > 10])
    _, hs = held(STEERING, prob, ["julia"])
    h = only(hs)
    @test h.resolved == 10
    @test Admission(:prerelease, "Stats", [12, 11]) in h.chain
    @test improved(h) == 12
    @test occursin("allow prereleases of Stats → allows: julia 12",
                   sprint(show, MIME("text/plain"), h))
end

@testset "holdback: when there is no remedy, say so" begin
    # Not every holdback has a fix. Here nothing the *user* set is to blame:
    # App's own bound on julia is what keeps julia back, and no constraint of
    # yours can move it -- only a release of App could. Inventing a fix would
    # be worse than admitting there is none.
    data = Dict(
        "App"   => PkgData([1], Dict(1 => ["julia"]),
                           Dict(1 => Dict("julia" => [10]))),
        "julia" => PkgData([12, 11, 10], HND(), HNC()),
    )
    prob = Problem(["App", "julia"])
    _, hs = held(data, prob, ["julia"])
    h = only(hs)
    @test h.resolved == 10
    @test h.best == 12
    @test h.fix === nothing
    @test improved(h) === nothing
    report = sprint(show, MIME("text/plain"), h)
    @test occursin("⇒ nothing you can change would move it.", report)
    @test occursin("nothing you can change would move it", summarize(h))
    @test sprint(show, h) == "Holdback(julia 10 < 12, no fix)"
end

@testset "holdback: best means best admissible" begin
    # "Below its best version" is measured against the versions this problem
    # would accept, not against the raw universe. A prerelease the query does
    # not admit is not an option being missed, so nothing is reported -- and
    # advertising it would bury the bound that is the real story.
    pre = Pair{Symbol,Any}[:prerelease => (p, v) -> p == "julia" && v == 12]
    prob = Problem(["App", "julia"]; excludes = pre)
    sol, hs = held(STEERING, prob, ["julia"])
    @test sol["julia"] == 11
    @test isempty(hs)      # 11 *is* the best julia this problem can reach

    # ... and with a bound that holds it below that best, the holdback is
    # against the best admissible version, not the prerelease above it
    prob = Problem(["App", "julia"];
                   compat = Dict("Stats" => [10]), excludes = pre)
    sol, hs = held(STEERING, prob, ["julia"])
    @test sol["julia"] == 10
    h = only(hs)
    @test h.best == 11
    @test improved(h) == 11
    @test occursin("julia resolved to 10; 11 is available.",
                   sprint(show, MIME("text/plain"), h))

    # a package with no admissible version at all is an unsatisfiable problem's
    # business; a holdback probe reports nothing rather than inventing a best
    none = Pair{Symbol,Any}[:prerelease => (p, v) -> p == "Stats"]
    prob = Problem(["App", "julia"]; compat = Dict("Stats" => [10]))
    info = prepare_pkg_info(pkg_info(STEERING, prob), prob)
    sol = resolve_prepared(info, prob)
    @test isempty(holdbacks(info, Problem(["App", "julia"]; excludes = none),
                            sol, ["Stats"]))
end

@testset "holdback: a truncated report says so" begin
    # Probing a package costs solves, so the budget is real -- and a report
    # that quietly dropped the rest would be misreporting. Ten packages, each
    # held back by its own compat bound, against a budget of three.
    data = Dict{String,PkgData{String,Int,Vector{Int},Vector{Int},HND,HNC}}()
    for i = 1:10
        data["P$i"] = PkgData([3, 2, 1], HND(), HNC())
    end
    reqs = ["P$i" for i = 1:10]
    prob = Problem(reqs; compat = Dict("P$i" => [2, 1] for i = 1:10))
    info = prepare_pkg_info(pkg_info(data, prob), prob)
    sol = resolve_prepared(info, prob)
    @test all(sol["P$i"] == 2 for i = 1:10)   # every one held below 3

    full = holdbacks(info, prob, sol; max_probes = 20)
    @test length(full) == 10
    @test full.unexamined == 0
    @test isempty(sprint(Resolver.render_unexamined, full))

    capped = holdbacks(info, prob, sol; max_probes = 3)
    @test length(capped) == 3
    @test capped.unexamined == 7
    # the ones it did examine are explained in full, as always
    @test all(h.fix !== nothing && improved(h) == 3 for h in capped)
    report = sprint(show, MIME("text/plain"), capped)
    @test occursin("(7 more packages resolved below their best version, not examined)",
                   report)
    # ... and the count is of candidates, not of confirmed holdbacks
    @test count("resolved to 2; 3 is available.", report) == 3

    # singular reads properly
    one = holdbacks(info, prob, sol; max_probes = 9)
    @test one.unexamined == 1
    @test occursin("(1 more package resolved below their best version, not examined)",
                   sprint(show, MIME("text/plain"), one))
end

@testset "holdback: bookkeeping" begin
    prob = Problem(["App", "julia"]; compat = Dict("Stats" => [10]))
    info = prepare_pkg_info(pkg_info(STEERING, prob), prob)
    sol = resolve_prepared(info, prob)

    # explaining does not disturb the universe it was handed, and repeats
    before = deepcopy(info)
    hs1 = holdbacks(info, prob, sol, ["julia"])
    @test info == before
    @test summarize(only(hs1)) == summarize(only(holdbacks(info, prob, sol, ["julia"])))
    # ... and the instance it builds is gone again, so a later resolve is
    # unaffected
    @test resolve_prepared(info, prob) == sol

    # every held-back package, and the probe budget that bounds it
    _, all = held(STEERING, prob)
    @test all isa Holdbacks{String,Int}
    @test sort!([h.pkg for h in all]) == ["Stats", "julia"]
    @test all.unexamined == 0
    _, capped = held(STEERING, prob; max_probes = 1)
    @test length(capped) == 1
    @test capped.unexamined == 1

    # asking about a package that is not in the solution, or not held back,
    # costs nothing and yields nothing
    @test isempty(holdbacks(info, prob, sol, ["App"]))
    @test isempty(holdbacks(info, prob, sol, ["Nonesuch"]))
end
