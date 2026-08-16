# Diagnosing an unsatisfiable resolve (src/Diagnostics.jl).
#
# Two oracles, neither of which asks the diagnosis anything:
#
#   * for the conflicts, the instance itself. A chain of facts claims that
#     assuming exactly those facts is unsatisfiable and that every one of them
#     is needed for that, and both halves are put to the solver directly.
#   * for the fixes, a resolve from the artifact. A fix claims that carrying
#     out its actions yields its solution, so the actions are turned into a
#     withdrawal here — independently of how the diagnosis turns them into one
#     — and the relaxed problem is prepared and resolved from scratch. That is
#     the path production dropped: the diagnosis answers on the universe the
#     failed resolve was run against, and this checks the two agree, which is
#     Theorem C where a report can see it. What the two share is `relax`'s
#     withdrawal of demands from a `Problem`, which is not this file's to check
#     (test/problem.jl and test/relaxation_stable.jl do).
#
# What the report says is checked separately, against hand-built cases, since
# there is no oracle for English.

using Resolver: Problem, PkgData, PkgInfo, SAT, Diagnosis, pkg_info, relax,
    prepare_pkg_info, finalize, is_satisfiable, installed_lit, forbidden_lit,
    with_classes_relaxed, class_exclusions, exclusion_kinds, nclasses,
    sat_assume_var
using Resolver.Diagnostics: Diagnostics, Conflict, Fix, Action, Fact,
    Requirement, Availability, unavailable, action_phrase

using Pkg.Versions: VersionSpec

# `random_problem` and `SHAPES` (relaxation.jl) are the constraint generator,
# and `verdicts` (relaxation.jl) the instance-comparison battery; that file
# runs before this one.

## helpers

# the instance a resolve of `prob` fails on, the universe behind it, and the
# artifact that universe was prepared from
function failed_instance(data, prob; order = nothing)
    info = pkg_info(data, prob)
    univ = prepare_pkg_info(info, prob; order)
    return SAT(univ), univ, info
end

# the literals a fact claims: the requirement it names, or every class this
# query emptied of the package it names. A fact is one package's worth of the
# story, so this is what taking the fact away would give back
fact_literals(sat::SAT{P}, f::Requirement{P}) where {P} =
    haskey(sat.vars, f.pkg) ? Int[installed_lit(sat, f.pkg)] : Int[]

fact_literals(sat::SAT{P}, f::Availability{P}) where {P} =
    haskey(sat.vars, f.pkg) ? Int[forbidden_lit(sat, f.pkg, c)
        for c in eachindex(sat.reps[f.pkg]) if iszero(sat.reps[f.pkg][c])] : Int[]

conflict_literals(sat::SAT, c::Conflict) =
    reduce(vcat, (fact_literals(sat, f) for f in c.chain); init = Int[])

# is the instance satisfiable assuming exactly `lits`?
function sat_assuming(sat::SAT, lits)
    for l in lits
        sat_assume_var(sat, l)
    end
    return is_satisfiable(sat)
end

# the withdrawal a set of actions asks for, read off here rather than taken
# from the diagnosis: `:drop` names a requirement to stop requiring, every
# other kind names a constraint of the problem and a package to lift it for
function withdrawal(actions::Vector{Action{P}}) where {P}
    drop_reqs = P[a.pkg for a in actions if a.kind === :drop]
    drop_constraints = Dict{Symbol,Set{P}}()
    for a in actions
        a.kind === :drop && continue
        push!(get!(Set{P}, drop_constraints, a.kind), a.pkg)
    end
    return drop_reqs, drop_constraints
end

# what carrying `actions` out gets you, the slow way: the relaxed problem
# prepared and resolved from the artifact, filter and all
fix_resolve(info, prob::Problem{P}, actions::Vector{Action{P}};
            order = nothing, by = identity) where {P} =
    resolve(info, relax(prob, withdrawal(actions)...); order, by, diagnose = false)

# every property a `Diagnosis` of `prob` claims, checked against the instance
# it was drawn from and against fresh resolves of the artifact. Returns the
# diagnosis, or `nothing` when the problem is satisfiable after all
function check_diagnosis(data, prob::Problem{P}; order = nothing,
                         by = identity) where {P}
    d = resolve(data, prob; order, by)
    d isa Diagnosis || return nothing
    sat, univ, info = failed_instance(data, prob; order)
    try
        ## the instance is left as it was found
        #
        # `diagnose` pops the deactivation frame, pushes a frame of definitional
        # clauses inside it, and resolves a relaxation per fix — so run it here
        # too and check that all of that came back. The checks below, and
        # anything else that reuses an instance, are entitled to find it as they
        # left it. The battery is one solve per class, so on a large universe
        # the verdict on the failed query is what is affordable
        battery = sum(nclasses(i) for i in values(sat.info); init = 0) ≤ 64
        before = battery ? verdicts(sat) : Bool[]
        deact = copy(sat.deact)
        @test Diagnostics.diagnose(sat, prob, univ; by, order) isa Diagnosis
        @test sat.deact == deact
        @test !issatisfiable(sat, prob.reqs)
        battery && @test verdicts(sat) == before

        ## conflicts

        # independent: no requirement is in two of them, and each names only
        # requirements
        seen = Set{P}()
        for c in d.conflicts
            @test !isempty(c.reqs)
            @test c.reqs ⊆ prob.reqs
            @test isempty(intersect(seen, c.reqs))
            union!(seen, c.reqs)
            # the chain names the conflict's requirements, then packages
            @test [f.pkg for f in c.chain if f isa Requirement] == c.reqs
            @test allunique(f.pkg for f in c.chain if f isa Availability)
        end

        with_classes_relaxed(sat) do
            for c in d.conflicts
                lits = conflict_literals(sat, c)
                if isempty(lits)
                    # a requirement the universe holds no version of: there is
                    # nothing to ask the instance, and the chain says as much
                    @test all(f -> !haskey(sat.vars, f.pkg), c.chain)
                    continue
                end
                # the chain is a conflict: assuming just these facts fails
                @test !sat_assuming(sat, lits)
                # ... and every fact in it is load-bearing: take one away and
                # what is left of the conflict is satisfiable
                for f in c.chain
                    drop = Set(fact_literals(sat, f))
                    @test sat_assuming(sat, filter(∉(drop), lits))
                end
            end
        end

        # availability facts say what the problem says about the versions, class
        # for class
        for c in d.conflicts, f in c.chain
            f isa Availability || continue
            if haskey(sat.info, f.pkg)
                info_p = sat.info[f.pkg]
                @test f.members == info_p.versions
                @test f.excluded ==
                    [exclusion_kinds(prob, f.pkg, v) for v in f.members]
                for k = 1:nclasses(info_p)
                    ex = class_exclusions(prob, f.pkg, info_p, k)
                    @test [f.members[i] for i in info_p.members[k]] == first.(ex)
                    @test [f.excluded[i] for i in info_p.members[k]] == last.(ex)
                end
                # a package with a version left over is not unavailable
                @test unavailable(f) == all(iszero, sat.reps[f.pkg])
            else
                # nothing of the package is in the universe to talk about
                @test isempty(f.members) && isempty(f.excluded)
                @test unavailable(f)
            end
        end

        ## fixes

        for fix in d.fixes
            @test !isempty(fix.actions)
            @test allunique(fix.actions)
            # the central property: carrying the actions out resolves, and to
            # exactly the solution the fix shows — the diagnosis answered the
            # relaxation on the universe filtered for the query, this answers it
            # on a universe filtered for the relaxation itself
            @test fix_resolve(info, prob, fix.actions; order, by) == fix.solution
        end
        # minimal and distinct: no fix asks for a strict superset of another's
        # actions, and no two ask for the same things
        sets = [Set(fix.actions) for fix in d.fixes]
        @test allunique(sets)
        for a in sets, b in sets
            @test !(b ⊊ a)
        end
    finally
        finalize(sat)
    end
    return d
end

## hand-built instances

const DEPS_NONE = Dict{Symbol,Vector{Symbol}}()
const COMP_NONE = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()

# A needs C@v1, B needs C@v2, and E and F disagree about G the same way: two
# conflicts that share nothing
const two_conflicts = Dict(
    :A => PkgData([:v1], Dict(:v1 => [:C]), Dict(:v1 => Dict(:C => [:v1]))),
    :B => PkgData([:v1], Dict(:v1 => [:C]), Dict(:v1 => Dict(:C => [:v2]))),
    :C => PkgData([:v2, :v1], DEPS_NONE, COMP_NONE),
    :E => PkgData([:v1], Dict(:v1 => [:G]), Dict(:v1 => Dict(:G => [:v1]))),
    :F => PkgData([:v1], Dict(:v1 => [:G]), Dict(:v1 => Dict(:G => [:v2]))),
    :G => PkgData([:v2, :v1], DEPS_NONE, COMP_NONE),
)

# one requirement, one dependency, and room for a bound to take the dependency
# away. Nothing in this registry tells :B's versions apart, so they are one
# class and it takes every one of them to empty it
const needs_dep = Dict(
    :A => PkgData([:v1], Dict(:v1 => [:B]), COMP_NONE),
    :B => PkgData([:w3, :w2, :w1], DEPS_NONE, COMP_NONE),
)

@testset "diagnosis: independent conflicts come out separately" begin
    d = check_diagnosis(two_conflicts, Problem([:A, :B, :E, :F]))
    @test d isa Diagnosis
    @test length(d.conflicts) == 2
    @test Set(Set(c.reqs) for c in d.conflicts) ==
        Set([Set([:A, :B]), Set([:E, :F])])
    # repairing needs one requirement out of each conflict: 2 × 2 ways
    @test length(d.fixes) == 4
    @test Set(Set(a.pkg for a in fix.actions) for fix in d.fixes) ==
        Set(Set([x, y]) for x in (:A, :B), y in (:E, :F))
    @test all(a.kind === :drop for fix in d.fixes for a in fix.actions)
    # and the requirements the conflict does not need stay out of it
    d = check_diagnosis(two_conflicts, Problem([:A, :B, :C]))
    @test length(d.conflicts) == 1
    @test d.conflicts[1].reqs == [:A, :B]
end

@testset "diagnosis: one conflict does not split" begin
    # three requirements that fail only together: each pair agrees on a version
    # of :D, all three on none
    data = Dict(
        :A => PkgData([:v1], Dict(:v1 => [:D]), Dict(:v1 => Dict(:D => [:v1, :v2]))),
        :B => PkgData([:v1], Dict(:v1 => [:D]), Dict(:v1 => Dict(:D => [:v2, :v3]))),
        :C => PkgData([:v1], Dict(:v1 => [:D]), Dict(:v1 => Dict(:D => [:v3, :v1]))),
        :D => PkgData([:v3, :v2, :v1], DEPS_NONE, COMP_NONE),
    )
    d = check_diagnosis(data, Problem([:A, :B, :C]))
    @test length(d.conflicts) == 1
    @test d.conflicts[1].reqs == [:A, :B, :C]
    # every pair of them is fine, which is what makes it one conflict
    for pair in ([:A, :B], [:A, :C], [:B, :C])
        @test resolve(data, Problem(pair); diagnose = false) !== nothing
    end
end

@testset "diagnosis: the story names the constraint" begin
    prob = Problem([:A]; compat = Dict(:B => [:w1]), pin = Dict(:B => :w2))
    d = check_diagnosis(needs_dep, prob)
    @test length(d.conflicts) == 1
    c = only(d.conflicts)
    @test c.reqs == [:A]
    @test c.chain == Fact[Requirement(:A),
        Availability{Symbol,Symbol}(:B, [:w3, :w2, :w1],
            [[:compat, :pin], [:compat], [:pin]])]
    @test unavailable(c.chain[2])
    # the cheapest version to give back costs one kind, so that is the fix —
    # and not requiring :A is the other
    @test [fix.actions for fix in d.fixes] ==
        [[Action(:compat, :B)], [Action(:drop, :A)]]
    @test d.fixes[1].solution == Dict(:A => :v1, :B => :w2)
    @test isempty(d.fixes[2].solution)
end

@testset "diagnosis: a package with a version left over" begin
    # :R@r1 needs :P and rules out :P@p1, and the bound takes :P@p2 away — so
    # :P is left with a version, just not one that works here
    data = Dict(
        :R => PkgData([:r1], Dict(:r1 => [:P]), Dict(:r1 => Dict(:P => [:p2]))),
        :P => PkgData([:p2, :p1], DEPS_NONE, COMP_NONE),
    )
    d = check_diagnosis(data, Problem([:R]; compat = Dict(:P => [:p1])))
    f = d.conflicts[1].chain[2]::Availability
    @test f.pkg == :P
    @test f.members == [:p2, :p1]
    @test f.excluded == [[:compat], Symbol[]]
    @test !unavailable(f)
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict, 2 fixes:

        Conflict 1: R cannot be satisfied.
          • you require R
          • P: p2 is excluded by your compat

        Verified fixes:
          1. relax your compat on P
             → allows: P p2, R r1
          2. drop requirement R
        """
end

@testset "diagnosis: any kind is a constraint like any other" begin
    prob = Problem([:A];
        prerelease = (p, v) -> p === :B && v !== :w1,
        yanked     = (p, v) -> p === :B && v === :w1)
    d = check_diagnosis(needs_dep, prob)
    f = d.conflicts[1].chain[2]::Availability
    @test f.excluded == [[:prerelease], [:prerelease], [:yanked]]
    # one kind is enough to give the package back, and the best version it
    # would give back costs the prerelease one
    @test [fix.actions for fix in d.fixes] ==
        [[Action(:prerelease, :B)], [Action(:drop, :A)]]
    @test d.fixes[1].solution == Dict(:A => :v1, :B => :w3)
    @test occursin("allow prerelease versions of B",
        sprint(show, MIME("text/plain"), d))
end

@testset "diagnosis: a requirement with nothing to install" begin
    data = Dict(
        :A => PkgData(Symbol[], DEPS_NONE, COMP_NONE),
        :B => PkgData([:v1], DEPS_NONE, COMP_NONE),
    )
    d = check_diagnosis(data, Problem([:A, :B]))
    @test length(d.conflicts) == 1
    @test d.conflicts[1].reqs == [:A]
    @test d.conflicts[1].chain == Fact[Requirement(:A),
        Availability{Symbol,Symbol}(:A, Symbol[], Vector{Symbol}[])]
    # nothing this query does is what took it away, so dropping it is all that
    # could help — and :B still resolves once it is gone
    @test [fix.actions for fix in d.fixes] == [[Action(:drop, :A)]]
    @test only(d.fixes).solution == Dict(:B => :v1)
    @test occursin("no version of A is available",
        sprint(show, MIME("text/plain"), d))
end

@testset "diagnosis: the instance is left as it was found" begin
    # Diagnosing takes the instance apart and puts it back: the deactivation
    # frame comes off, a frame of definitional clauses goes on inside it, both
    # are undone, and then a relaxation is resolved on the instance per fix. An
    # instance is reusable, so what has to be true afterwards is that it answers
    # every question exactly as it did — the emptied classes still emptied,
    # every verdict what it was.
    #
    # What the instance answers is the whole of the contract, and it is the
    # only thing worth asserting: the solver's clause and variable counts both
    # grow, since popping a frame satisfies its clauses rather than removing
    # them and the definitional variables outlive the clauses that defined
    # them. Neither moves an answer, and a bare lift with no diagnosis at all
    # grows the clause count exactly the same way.
    unavailable_dep = Problem([:A];
        compat = Dict(:B => [:w1]), pin = Dict(:B => :w2))
    one_left = Problem([:R]; compat = Dict(:P => [:p1]))
    left_over = Dict(
        :R => PkgData([:r1], Dict(:r1 => [:P]), Dict(:r1 => Dict(:P => [:p2]))),
        :P => PkgData([:p2, :p1], DEPS_NONE, COMP_NONE),
    )
    cases = [
        # a query that emptied classes of one package ...
        needs_dep => unavailable_dep,
        left_over => one_left,
        # ... of several ...
        needs_dep => Problem([:A]; compat = Dict(:A => Symbol[], :B => Symbol[])),
        # ... and one that emptied none at all, so there is no frame to lift
        two_conflicts => Problem([:A, :B, :E, :F]),
    ]
    for (data, prob) in cases
        sat, univ, _ = failed_instance(data, prob)
        try
            before = verdicts(sat)
            deact = copy(sat.deact)
            # diagnosing twice over, since restoring has to be repeatable
            for _ = 1:2
                @test Diagnostics.diagnose(sat, prob, univ) isa Diagnosis
                @test verdicts(sat) == before
                @test sat.deact == deact
                @test !issatisfiable(sat, prob.reqs)
            end
            # and the frame is genuinely back rather than never having been
            # there: re-imposing it by assumption says what it says
            @test verdicts(sat, sat.deact) == before
        finally
            finalize(sat)
        end
    end
end

## the API

@testset "diagnosis: what resolve returns" begin
    prob = Problem([:A]; compat = Dict(:B => Symbol[]))
    for src in (needs_dep, pkg_info(needs_dep, prob))
        @test resolve(src, Problem([:A])) isa Dict{Symbol,Symbol}
        @test resolve(src, prob) isa Diagnosis{Symbol,Symbol}
        @test resolve(src, prob; diagnose = false) === nothing
        # ... and through the bare-requirements form, which is the
        # unconstrained problem and so diagnoses the same way
        @test resolve(src, [:A]) isa Dict{Symbol,Symbol}
    end
    # `issatisfiable` is untouched by any of it
    @test !issatisfiable(needs_dep, prob)
    @test issatisfiable(needs_dep, Problem([:A]))
    # a caller-supplied T1 artifact comes through a diagnosed resolve intact,
    # as it does through any other
    info = pkg_info(needs_dep, prob)
    before = deepcopy(info)
    @test resolve(info, prob) isa Diagnosis
    @test info == before
    # ... and the data-dict entry point, whose artifact the resolve is allowed
    # to consume, diagnoses just the same
    @test resolve(copy(needs_dep), prob) isa Diagnosis
end

## rendering

@testset "diagnosis: the report" begin
    d = resolve(two_conflicts, Problem([:A, :B, :E, :F]))
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 2 conflicts, 4 fixes:

        Conflict 1: A and B cannot be satisfied together.
          • you require A
          • you require B

        Conflict 2: E and F cannot be satisfied together.
          • you require E
          • you require F

        Verified fixes:
          1. drop requirement B and drop requirement F
             → allows: A v1, E v1
          2. drop requirement B and drop requirement E
             → allows: A v1, F v1
          3. drop requirement A and drop requirement F
             → allows: B v1, E v1
          4. drop requirement A and drop requirement E
             → allows: B v1, F v1
        """
    # the one-line summary
    @test sprint(show, d) == "Diagnosis: 2 conflicts, 4 fixes"

    # every kind that excludes a version is named, and a long line is filled
    d = resolve(needs_dep,
        Problem([:A]; compat = Dict(:B => [:w1]), pin = Dict(:B => :w2)))
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict, 2 fixes:

        Conflict 1: A cannot be satisfied.
          • you require A
          • no version of B is available: w3 is excluded by your compat and your pin,
            w2 by your compat, w1 by your pin

        Verified fixes:
          1. relax your compat on B
             → allows: A v1, B w2
          2. drop requirement A
        """
end

@testset "diagnosis: the report says nothing about how it was found" begin
    # classes, literals, assumptions, cores and the solver are how the answer
    # was found; what it says is about packages, versions and constraints
    cases = [
        two_conflicts => Problem([:A, :B, :E, :F]),
        needs_dep => Problem([:A]; compat = Dict(:B => Symbol[])),
        needs_dep => Problem([:A]; pin = Dict(:B => :w9)),
        needs_dep => Problem([:A]; compat = Dict(:B => Symbol[]),
                                   pin = Dict(:A => :v9)),
    ]
    for (data, prob) in cases
        d = resolve(data, prob)
        @test d isa Diagnosis
        report = sprint(show, MIME("text/plain"), d)
        for word in ("class", "literal", "assum", "MUS", "MCS", "core",
                     "solver", "SAT", "clause", "variable", "deactivat")
            @test !occursin(word, report)
        end
    end
    # an action reads as something the user could carry out, whatever the kind
    # of the constraint it lifts is called
    @test action_phrase(Action(:drop, :A)) == "drop requirement A"
    @test action_phrase(Action(:compat, :A)) == "relax your compat on A"
    @test action_phrase(Action(:pin, :A)) == "unpin A"
    @test action_phrase(Action(:prerelease, :A)) ==
        "allow prerelease versions of A"
end

## sweeps

@testset "diagnosis: verified fixes over generated data" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    diagnoses = fixes = relaxations = 0
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:10
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            for shape in SHAPES, by in (identity, p -> -p)
                prob = random_problem(reqs, m, n, shape)
                diag = check_diagnosis(data, prob; by)
                diag === nothing && continue
                diagnoses += 1
                fixes += length(diag.fixes)
                relaxations += count(a.kind !== :drop
                    for fix in diag.fixes for a in fix.actions)
            end
        end
    end
    # the sweep really did diagnose, the diagnoses really did offer fixes, and
    # the fixes really did include relaxing a constraint rather than only
    # dropping requirements
    @test diagnoses > 0
    @test fixes > diagnoses
    @test relaxations > 0
end

@testset "diagnosis: verified fixes under a version ordering" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    diagnosed = 0
    up = p -> (u, v) -> u > v # prefer the lowest version
    for (m, n) in ((2, 3), (3, 2), (3, 3))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:8
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            prob = random_problem(reqs, m, n, :all)
            for order in (nothing, up)
                diagnosed += check_diagnosis(data, prob; order) !== nothing
            end
        end
    end
    @test diagnosed > 0
end

@testset "diagnosis: registry-scale" begin
    rp = registry.provider()
    # a bound that leaves DataFrames' table printer with nothing, on a
    # DataFrames new enough to need one
    prob = Problem(["DataFrames"];
        compat = Dict("DataFrames" => VersionSpec("1"),
                      "PrettyTables" => VersionSpec("99")))
    d = check_diagnosis(rp, prob)
    @test d isa Diagnosis{String,VersionNumber}
    c = only(d.conflicts)
    @test c.reqs == ["DataFrames"]
    # the story is the requirement, the bound that forces a modern DataFrames,
    # and the bound that leaves it without a table printer
    @test any(f isa Availability && f.pkg == "PrettyTables" && unavailable(f)
              for f in c.chain)
    # relaxing either bound is a fix, and so is not requiring DataFrames
    sets = Set(Set(fix.actions) for fix in d.fixes)
    @test Set([Action(:compat, "PrettyTables")]) ∈ sets
    @test Set([Action(:drop, "DataFrames")]) ∈ sets
    for fix in d.fixes
        fix.actions == [Action(:compat, "PrettyTables")] || continue
        @test fix.solution["DataFrames"] ∈ VersionSpec("1")
        @test haskey(fix.solution, "PrettyTables")
    end
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("no version of PrettyTables is available", report)
    @test occursin("relax your compat on PrettyTables", report)
end
