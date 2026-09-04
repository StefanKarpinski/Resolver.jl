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
# A report hands out one menu per conflict and lets the reader choose from each
# independently, so the claim it rests on is that *every* combination of choices
# repairs the query. The menus are small, so every check_diagnosis checks every
# combination, from scratch. What the combinations leave out (`others`) is
# checked against a third oracle, in the cases small enough for it: every subset
# of a pool of actions, resolved from the artifact, reduced to the minimal ones.
#
# Half of what `others` says is that no repair costs more than the ones on the
# menus, which the diagnosis decides with a single solve against a bounded
# enumeration. That gets a fourth oracle: the whole family of repairs,
# enumerated with no bound at all, judged by its sizes.
#
# What the report says is checked separately, against hand-built cases, since
# there is no oracle for English.

using Resolver: Problem, PkgData, PkgInfo, SAT, Diagnosis, pkg_info, relax,
    prepare_pkg_info, finalize, sat_solve, installed_lit, forbidden_lit,
    with_classes_relaxed, class_exclusions, exclusion_kinds, nclasses,
    sat_assume_var
using Resolver.Diagnostics: Diagnostics, Conflict, Fix, Action, Line,
    clause_versions, clauses_satisfiable, clause_of, project, action_phrase
using Resolver.Clauses: Clauses, Clause, packages, isbottom, clause_phrase,
    literal, resolve_on, subsumes
using Resolver.UnsatCores: sat_mcses

@isdefined(ProofCheck) || include(joinpath(@__DIR__, "proof_check.jl"))
using .ProofCheck
using .ProofCheck: printed_lines

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

# is the instance satisfiable assuming exactly `lits`?
function sat_assuming(sat::SAT, lits)
    for l in lits
        sat_assume_var(sat, l)
    end
    return sat_solve(sat)
end



# the versions of `p` a fact speaks for, as indices into what `p` offers
fact_indices(sat::SAT{P}, p::P, versions) where {P} =
    Int[findfirst(==(v), sat.info[p].versions) for v in versions]




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

# one fix out of every conflict's menu, as the actions to carry out together
# (over copies, so that a one-conflict choice does not hand back the fix's own
# vector for `unique!` to work on)
combined(choice) = unique!(reduce(vcat, (copy(fix.actions) for fix in choice)))

# every way of choosing one fix per conflict, as a set of action sets — the
# fixes of the whole query the report offers
fix_combinations(d::Diagnosis) =
    Set(Set(combined(choice))
        for choice in Iterators.product((c.fixes for c in d.conflicts)...))

# every minimal repair drawn from `pool`, the slow way: each subset of it
# prepared and resolved from the artifact, and the ones that resolve reduced to
# the minimal ones. The oracle for what a report's combinations leave out
function minimal_repairs(info, prob::Problem{P}, pool::Vector{Action{P}};
                         order = nothing, by = identity) where {P}
    repairs = Set{Set{Action{P}}}()
    for bits = 0:(1 << length(pool)) - 1
        actions = Action{P}[pool[i] for i in eachindex(pool)
                            if !iszero(bits & (1 << (i-1)))]
        fix_resolve(info, prob, actions; order, by) === nothing && continue
        push!(repairs, Set(actions))
    end
    return Set(r for r in repairs if !any(s -> s ⊊ r, repairs))
end

# How many repairs to enumerate for the oracle below before giving up on it.
const REPAIR_LIMIT = 512

# Every repair of `prob` there is, as sets of the literals a diagnosis assumes,
# or `nothing` when there are too many to enumerate. The oracle for what a
# report says it leaves out: the diagnosis decides whether anything costs more
# than its menus do with a single solve, and this decides it by enumerating the
# whole family and looking at the sizes.
#
# It asks the same questions of the same instance `diagnose` does, so it goes
# through the same helpers — and, like `diagnose`, leaves the instance as it
# found it.
function all_repairs(sat::SAT{P}, prob::Problem{P}) where {P}
    vm = Diagnostics.VarMap(sat)
    req_lits = Int[installed_lit(sat, p) for p in unique(prob.reqs)
                   if haskey(vm, p)]
    repairs = with_classes_relaxed(sat) do
        Diagnostics.with_emptied_packages(sat, vm) do pkg_lits, _
            sat_mcses(sat, [req_lits; pkg_lits]; limit = REPAIR_LIMIT)
        end
    end
    return length(repairs) < REPAIR_LIMIT ? repairs : nothing
end

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

        # each names requirements of the query
        for c in d.conflicts
            @test !isempty(c.reqs)
            @test c.reqs ⊆ prob.reqs
        end

        # A proof, checked as two things and no more: every statement on the
        # page is true of the universe this query left, and they cannot hold
        # together. A set of true statements that contradict is the whole of
        # what proving unsatisfiability is.
        #
        # There is no third question about how one line follows from another,
        # because none of them does: the lines are flat, each derived from the
        # registry rather than from its neighbours. What used to be checked
        # besides -- that a line spoke in a direction the registry licensed,
        # that a bound was stated where one was claimed -- is not a separate
        # question either. A clause has no direction to get wrong, and its
        # bound is the literal.
        for c in d.conflicts
            problems = proof_problems(sat, prob, c)
            @test isempty(problems) || (@show problems; false)
        end

        # ... and within one proof no line says less than another beside it.
        # Across proofs it may: each is an argument on its own and has to stand
        # without borrowing a step, so a line two of them both need is written
        # out in both. Dropping a line need not restore satisfiability either,
        # for the same reason -- another proof is still standing.
        for c in d.conflicts
            for n in unique!(Int[l.proof for l in c.lines])
                cs = Clause{P}[l.clause for l in c.lines if l.proof == n]
                for a in cs, b in cs
                    a == b && continue
                    @test !subsumes(a, b) ||
                        (@show ("a line is covered", n, a, b); false)
                end
            end
        end

        # ... and it accounts for the menu under it. A fix names packages to
        # act on, and a report that offered one it never spoke of would be
        # talking past itself: the reader is told to change something the
        # argument never mentioned
        for c in d.conflicts
            named = Set{P}(p for l in c.lines for p in packages(l.clause))
            for fix in c.fixes, a in fix.actions
                isempty(c.lines) && continue   # a package with no versions
                @test a.pkg in named
            end
        end

        # the versions a conflict carries are the ones its statements are about
        for c in d.conflicts
            @test Set(keys(c.versions)) ==
                  Set(p for l in c.lines for p in packages(l.clause))
        end

        ## fixes

        for (i, c) in enumerate(d.conflicts)
            for fix in c.fixes
                @test !isempty(fix.actions)
                @test allunique(fix.actions)
                # a fix shows what it gets you when the other conflicts are
                # fixed the first way their menus offer, so that is what is
                # resolved here — the diagnosis answered the relaxation on the
                # universe filtered for the query, this answers it on a
                # universe filtered for the relaxation itself
                actions = combined([j == i ? fix : first(d.conflicts[j].fixes)
                                    for j in eachindex(d.conflicts)])
                @test fix_resolve(info, prob, actions; order, by) == fix.solution
            end
            # minimal and distinct: within one menu, no fix asks for a strict
            # superset of another's actions, and no two ask for the same things
            sets = [Set(fix.actions) for fix in c.fixes]
            @test allunique(sets)
            for a in sets, b in sets
                @test !(b ⊊ a)
            end
        end

        ## every combination is a fix
        #
        # The claim the whole report rests on: the menus are independent, so
        # any one fix from each of them repairs the query. Checked from
        # scratch, exhaustively — the menus are small enough that they can be
        menus = [c.fixes for c in d.conflicts]
        if !isempty(menus) && prod(length, menus) ≤ 64
            for choice in Iterators.product(menus...)
                @test fix_resolve(info, prob, combined(choice);
                    order, by) !== nothing
            end
        end

        ## what the menus leave out
        #
        # Half of what `others` says is that no repair of the query costs more
        # than the ones the menus reach. The diagnosis decides that by asking
        # the instance for a repair holding none of the cheapest ones; this
        # decides it by enumerating every repair there is and comparing sizes.
        # `:some` claims nothing either way — it is what a family the menus do
        # not cover reports, whatever else is true of it
        repairs = all_repairs(sat, prob)
        if repairs !== nothing && !isempty(repairs)
            larger = maximum(length, repairs) > minimum(length, repairs)
            d.others === :none   && @test !larger
            d.others === :larger && @test larger
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

# two requirements, one dependency between them: giving :C back fixes both at
# once, and dropping both requirements is the only other way — a repair that
# gives up strictly more, so a report that offers the first has left it out
const shared_dep = Dict(
    :A => PkgData([:v1], Dict(:v1 => [:C]), COMP_NONE),
    :B => PkgData([:v1], Dict(:v1 => [:C]), COMP_NONE),
    :C => PkgData([:c1], DEPS_NONE, COMP_NONE),
)

# two requirements that each need :C at its newer version, and a bound of their
# own holding each of them at the version that does. One bound on :C is what
# breaks both, and relaxing it is the one cheapest fix — so the conflict has to
# name both requirements, not whichever of them the search reached first
const shared_bound = Dict(
    :A => PkgData([:a2, :a1], Dict(:a2 => [:C]), Dict(:a2 => Dict(:C => [:c2]))),
    :B => PkgData([:b2, :b1], Dict(:b2 => [:C]), Dict(:b2 => Dict(:C => [:c2]))),
    :C => PkgData([:c2, :c1], DEPS_NONE, COMP_NONE),
)

# :A and :B both need :C, and :E and :F disagree about :G. With a bound that
# leaves :C nothing, settling the first conflict cheaply takes one action and
# settling it any other way takes two — so the smallest repairs are a product of
# two menus and there are repairs beyond them that cost more
const larger_repair = Dict(
    :A => PkgData([:v1], Dict(:v1 => [:C]), COMP_NONE),
    :B => PkgData([:v1], Dict(:v1 => [:C]), COMP_NONE),
    :C => PkgData([:c1], DEPS_NONE, COMP_NONE),
    :E => PkgData([:v1], Dict(:v1 => [:G]), Dict(:v1 => Dict(:G => [:v1]))),
    :F => PkgData([:v1], Dict(:v1 => [:G]), Dict(:v1 => Dict(:G => [:v2]))),
    :G => PkgData([:v2, :v1], DEPS_NONE, COMP_NONE),
)

# :A's versions disagree about :C: a3 wants c3, the older two want c2. A bound
# leaving only c1 breaks all of them, and for different reasons — so a story
# that gave one bound for the package would be giving one none of its versions
# declares
const varied_bound = Dict(
    :A => PkgData([:a3, :a2, :a1],
        Dict(:a3 => [:C], :a2 => [:C], :a1 => [:C]),
        Dict(:a3 => Dict(:C => [:c3]), :a2 => Dict(:C => [:c2]),
             :a1 => Dict(:C => [:c2]))),
    :C => PkgData([:c3, :c2, :c1], DEPS_NONE, COMP_NONE),
)

# :A needs :B needs :C, and the bound is on :C. Nothing relates :A to :C, so a
# story drawn only from the packages the chain names would have nothing to say
const two_hops = Dict(
    :A => PkgData([:a1], Dict(:a1 => [:B]), COMP_NONE),
    :B => PkgData([:b1], Dict(:b1 => [:C]), Dict(:b1 => Dict(:C => [:c2]))),
    :C => PkgData([:c2, :c1], DEPS_NONE, COMP_NONE),
)

# :P names :W in its compat and does not depend on it — a weak dependency,
# which is a bound with no edge behind it. :S is what pulls :W in
const weak_bound = Dict(
    :P => PkgData([:p1], DEPS_NONE, Dict(:p1 => Dict(:W => [:w1]))),
    :S => PkgData([:s1], Dict(:s1 => [:W]), COMP_NONE),
    :W => PkgData([:w2, :w1], DEPS_NONE, COMP_NONE),
)

# four requirements whose disagreements form a path — :C with :B, :B with :A,
# :A with :D — so the cheapest repairs are {A,B}, {A,C} and {B,D}, which is
# three of them and therefore not a product of anything
const conflict_path = Dict(
    :A => PkgData([:v1], Dict(:v1 => [:P, :Q]),
        Dict(:v1 => Dict(:P => [:p1], :Q => [:q1]))),
    :B => PkgData([:v1], Dict(:v1 => [:P, :R]),
        Dict(:v1 => Dict(:P => [:p2], :R => [:r1]))),
    :C => PkgData([:v1], Dict(:v1 => [:R]), Dict(:v1 => Dict(:R => [:r2]))),
    :D => PkgData([:v1], Dict(:v1 => [:Q]), Dict(:v1 => Dict(:Q => [:q2]))),
    :P => PkgData([:p2, :p1], DEPS_NONE, COMP_NONE),
    :Q => PkgData([:q2, :q1], DEPS_NONE, COMP_NONE),
    :R => PkgData([:r2, :r1], DEPS_NONE, COMP_NONE),
)

# Two conflicts at once, the second of them squeezed twice over: :B needs :U at
# a version the query has taken away, and :E needs :Z at one and :T at another,
# reaching :T through :M. Either of :E's squeezes explains it by itself, so the
# smallest explanation names one of them — and the walk from :E meets the other
# on its way there
const two_squeezes = Dict(
    :B => PkgData([:b1], Dict(:b1 => [:U]), Dict(:b1 => Dict(:U => [:u2]))),
    :E => PkgData([:e2, :e1], Dict(:e1 => [:M, :Z]),
        Dict(:e1 => Dict(:Z => [:z2]))),
    :M => PkgData([:m1], Dict(:m1 => [:T]), Dict(:m1 => Dict(:T => [:t2]))),
    :T => PkgData([:t2, :t1], DEPS_NONE, COMP_NONE),
    :U => PkgData([:u2, :u1], DEPS_NONE, COMP_NONE),
    :Z => PkgData([:z2, :z1], DEPS_NONE, COMP_NONE),
)

# :A needs :C at a version the query has taken away, and what is left of :C
# needs :A back. The story ends at :C, so :C is a package the chain speaks
# *about* before it ever speaks *for* it — and what the query left of it is a
# premise of the last line rather than of the first
const late_speaker = Dict(
    :A => PkgData([:a1], Dict(:a1 => [:C]), Dict(:a1 => Dict(:C => [:c2]))),
    :C => PkgData([:c2, :c1], Dict(:c1 => [:A]), COMP_NONE),
)

# :A's newest version would have done, and the query took it away: the chain
# says so, and the run that is left is what that sentence answers
const narrowed_run = Dict(
    :A => PkgData([:a3, :a2, :a1],
        Dict(:a3 => [:B], :a2 => [:B], :a1 => [:B]),
        Dict(:a3 => Dict(:B => [:b1]), :a2 => Dict(:B => [:b2]),
             :a1 => Dict(:B => [:b2]))),
    :B => PkgData([:b2, :b1], DEPS_NONE, COMP_NONE),
)

# every version of :A there is agrees about :B, and the query takes none of
# them away: there is nothing for a range over :A to be answering
const whole_run = Dict(
    :A => PkgData([:a2, :a1], Dict(:a2 => [:B], :a1 => [:B]),
        Dict(:a2 => Dict(:B => [:b2]), :a1 => Dict(:B => [:b2]))),
    :B => PkgData([:b2, :b1], DEPS_NONE, COMP_NONE),
)

@testset "diagnosis: independent conflicts come out separately" begin
    prob = Problem([:A, :B, :E, :F])
    d = check_diagnosis(two_conflicts, prob)
    @test d isa Diagnosis
    @test length(d.conflicts) == 2
    @test Set(Set(c.reqs) for c in d.conflicts) ==
        Set([Set([:A, :B]), Set([:E, :F])])
    # each conflict is settled by dropping one of its own requirements, and
    # says nothing about the other's: 2 and 2, not 2 × 2
    for c in d.conflicts
        @test Set(Set(a.pkg for a in fix.actions) for fix in c.fixes) ==
            Set(Set([p]) for p in c.reqs)
        @test all(a.kind === :drop for fix in c.fixes for a in fix.actions)
    end
    # ... and one from each is a repair, in all 2 × 2 ways. There is nothing
    # else: the combinations are exactly the minimal repairs, which is what
    # entitles the report to say nothing further
    pool = [Action(:drop, p) for p in (:A, :B, :E, :F)]
    @test fix_combinations(d) ==
        minimal_repairs(pkg_info(two_conflicts, prob), prob, pool)
    @test d.others === :none
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

@testset "a proof answers to its own claim, not to the union" begin
    # A chain is a union of proofs, one per thing the menu offers to undo, and
    # "the union is still unsatisfiable" is far too weak a test of any change
    # to it. Here two independent claims share a chain: `:A` needs `:B` and the
    # query leaves `:B` nothing; `:C` needs `:D` and the query leaves `:D`
    # nothing. Either alone makes the union unsatisfiable -- so a check against
    # the union would happily delete the whole of the other one, leaving a fix
    # the report can no longer account for.
    data = Dict(
        :A => PkgData([:a1], Dict(:a1 => [:B]), COMP_NONE),
        :B => PkgData([:b1], DEPS_NONE, COMP_NONE),
        :C => PkgData([:c1], Dict(:c1 => [:D]), COMP_NONE),
        :D => PkgData([:d1], DEPS_NONE, COMP_NONE),
    )
    prob = Problem([:A, :C]; compat = Dict(:B => Symbol[], :D => Symbol[]))
    sat, _, _ = failed_instance(data, prob)
    try
        # Each claim's proof is derived on its own, so the report has to
        # account for both: a statement about :A needing :B, and one about :C
        # needing :D. Derived from the union instead, either claim's lines
        # could be thrown away whole -- the union stays unsatisfiable on the
        # strength of the other -- and the menu would offer a fix the argument
        # never speaks of.
        d = resolve(data, prob)
        # two reasons, so two conflicts -- and each states its own, rather
        # than one of them riding on the other's unsatisfiability
        @test length(d.conflicts) == 2
        all_named = Set{Symbol}()
        for c in d.conflicts
            named = Set(p for l in c.lines for p in packages(l.clause))
            union!(all_named, named)
            # the fixes this conflict offers are the ones its own proof speaks of
            for fix in c.fixes, act in fix.actions
                @test act.pkg in named
            end
        end
        @test Set([:A, :B, :C, :D]) ⊆ all_named
    finally
        finalize(sat)
    end
end

@testset "diagnosis: the story names the constraint" begin
    prob = Problem([:A]; compat = Dict(:B => [:w1]), pin = Dict(:B => :w2))
    d = check_diagnosis(needs_dep, prob)
    @test length(d.conflicts) == 1
    c = only(d.conflicts)
    @test c.reqs == [:A]
    # the middle fact is the one the user did not write: :A needs :B, and no
    # bound of :A's is what the story turns on, since the query has left :B
    # with nothing whatever :A would have taken
    # the cheapest version to give back costs one kind, so that is the fix —
    # and not requiring :A is the other
    @test [fix.actions for fix in c.fixes] ==
        [[Action(:compat, :B)], [Action(:drop, :A)]]
    @test c.fixes[1].solution == Dict(:A => :v1, :B => :w2)
    @test isempty(c.fixes[2].solution)
    @test d.others === :none
end

@testset "diagnosis: a package with a version left over" begin
    # :R@r1 needs :P and rules out :P@p1, and the bound takes :P@p2 away — so
    # :P is left with a version, just not one that works here
    data = Dict(
        :R => PkgData([:r1], Dict(:r1 => [:P]), Dict(:r1 => Dict(:P => [:p2]))),
        :P => PkgData([:p2, :p1], DEPS_NONE, COMP_NONE),
    )
    d = check_diagnosis(data, Problem([:R]; compat = Dict(:P => [:p1])))
    # here there is a bound to state: :P has a version left, just not one
    # :R will take
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict:

        Conflict 1: R cannot be satisfied.
          • you require R
          • your compat leaves P p1
          • R requires P p2
          Fix it by any one of:
            1. relax your compat on P
               → allows: P p2, R r1
            2. drop requirement R
        """
end

@testset "diagnosis: a bound that differs across versions" begin
    # :C has a version left, and every version of :A wants a different one it
    # has not got. What each of them wants is a fact about that version, so
    # the story states it once per version that agrees
    d = check_diagnosis(varied_bound, Problem([:A]; compat = Dict(:C => [:c1])))
    c = only(d.conflicts)
    # :C is the only package this chain links :A to, so nothing it states can
    # tell :a3 from :a2. The cases are resolved on :A in one go -- intersect
    # what each speaks for, union what each forces -- and since between them
    # they speak for every version, what is left is one line about what any
    # of them forces
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict:

        Conflict 1: A cannot be satisfied.
          • you require A
          • your compat leaves C c1
          • A requires C ≥c2
          Fix it by any one of:
            1. relax your compat on C
               → allows: A a3, C c3
            2. drop requirement A
        """
end

@testset "diagnosis: a package's availability is a premise where it speaks" begin
    # :C is the package the query emptied and the package the story ends at,
    # and it also speaks: what is left of it needs :A back. So the chain says
    # what the query left of it where that becomes a premise — after the line
    # that asks for it and before the line that speaks for it. Hoisting every
    # availability the relations name to the front puts this one first, where
    # nothing has introduced :C yet and it reads as arriving from nowhere
    d = check_diagnosis(late_speaker, Problem([:A]; compat = Dict(:C => [:c1])))
    c = only(d.conflicts)
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict:

        Conflict 1: A cannot be satisfied.
          • you require A
          • your compat leaves C c1
          • A requires C c2
          Fix it by any one of:
            1. relax your compat on C
               → allows: A a1, C c2
            2. drop requirement A
        """
end

@testset "diagnosis: a subject that is the whole package is not a range" begin
    # Nothing has narrowed :A — the query takes none of its versions away, and
    # the chain says nothing about what it left of it — so naming the run in
    # full would read as a narrowing that never happened. The package is what
    # the line is about, and the package is what it says
    d = check_diagnosis(whole_run, Problem([:A]; compat = Dict(:B => [:b1])))
    c = only(d.conflicts)
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("requires B b2", report)
    @test !occursin("a1–a2", report)
    # ... where the chain does say what the query left of the package, the run
    # is what that sentence is answering, and it is named in full
    d2 = check_diagnosis(narrowed_run,
        Problem([:A]; compat = Dict(:A => [:a2, :a1], :B => [:b1])))
    report2 = sprint(show, MIME("text/plain"), d2)
    @test occursin("your compat leaves ≤a2", report2)
    @test occursin("A ≤a2 requires B b2", report2)
end

@testset "diagnosis: a story that spans more than one hop" begin
    # the query names :A and :C, and it is :B in between that explains them
    d = check_diagnosis(two_hops, Problem([:A]; compat = Dict(:C => [:c1])))
    c = only(d.conflicts)
    # :B is resolved away, so the line is about the two packages the query
    # named and says which package it reached them through
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("A requires C c2 (through B)", report)
    @test !occursin("B requires", report)
    # ... and :B is not something the query said anything about: the middle of
    # the story is the part only the registry knows
    # ... and :B is the middle of the story: the query said nothing whatever
    # about it, so it is neither required nor constrained -- the proof reaches
    # it through the registry alone
    @test :B ∉ c.reqs
    @test :B ∉ keys(c.excluded)
end

@testset "diagnosis: a bound with no dependency behind it" begin
    # :P's bound on :W is a weak dependency: a compat entry with no edge. The
    # universe records compatibility symmetrically, so which of the two
    # declared it is not recoverable — and saying ":P requires :W" would be
    # attributing a bound to a package that does not even depend on the other
    prob = Problem([:P, :S]; compat = Dict(:W => [:w2]))
    d = check_diagnosis(weak_bound, prob)
    c = only(d.conflicts)
    # each requirement's own consequence follows it: the chain is stored in
    # the order it argues in
    # the only dependency stated is the one the registry has
    # the only dependency stated is the one the registry has: :P's bound on
    # :W permits :W's absence, so nothing on the page says :P brings it in
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("S requires W", report)
    @test occursin("P constrains W w1", report)
    @test !occursin("P requires", report)
    # ... and the two lines do not leave :W nothing on their own -- what rules
    # out :w1 is the query -- so the report does not claim that they do
    @test !occursin("all of these", report)
    # ... and it really is a weak dependency: :P on its own installs no :W, so
    # nothing about it can be a dependency of :P
    sol = resolve(weak_bound, [:P])
    @test sol !== nothing && !haskey(sol, :W)
end

@testset "diagnosis: any kind is a constraint like any other" begin
    prob = Problem([:A];
        prerelease = (p, v) -> p === :B && v !== :w1,
        yanked     = (p, v) -> p === :B && v === :w1)
    d = check_diagnosis(needs_dep, prob)
    # one kind is enough to give the package back, and the best version it
    # would give back costs the prerelease one
    fixes = only(d.conflicts).fixes
    @test [fix.actions for fix in fixes] ==
        [[Action(:prerelease, :B)], [Action(:drop, :A)]]
    @test fixes[1].solution == Dict(:A => :v1, :B => :w3)
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
    # a package with no versions has a domain of one element, so nothing can
    # be said about it and there is no proof to print -- the heading and the
    # one fix are the whole of it
    @test isempty(d.conflicts[1].lines)
    # nothing this query does is what took it away, so dropping it is all that
    # could help — and :B still resolves once it is gone
    fixes = only(d.conflicts).fixes
    @test [fix.actions for fix in fixes] == [[Action(:drop, :A)]]
    @test only(fixes).solution == Dict(:B => :v1)
    @test d.others === :none
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("no version of A is available", report)
    # one thing to do reads as one thing to do, not as a menu of one
    @test occursin("The only fix: drop requirement A", report)
    @test !occursin("any one of", report)
end

@testset "diagnosis: a repair that gives up more is not on the menu" begin
    # :A and :B both need :C, and the bound leaves :C with nothing. Giving :C
    # back is the one cheapest repair; dropping both requirements is a repair
    # too, and gives up strictly more, so the report says there is more
    prob = Problem([:A, :B]; compat = Dict(:C => Symbol[]))
    d = check_diagnosis(shared_dep, prob)
    c = only(d.conflicts)
    @test [fix.actions for fix in c.fixes] == [[Action(:compat, :C)]]
    @test only(c.fixes).solution == Dict(:A => :v1, :B => :v1, :C => :c1)
    # the oracle: the minimal repairs are those two, and the report offers the
    # smaller of them and nothing else
    info = pkg_info(shared_dep, prob)
    pool = [Action(:compat, :C), Action(:drop, :A), Action(:drop, :B)]
    repairs = minimal_repairs(info, prob, pool)
    @test repairs == Set([Set([Action(:compat, :C)]),
                          Set([Action(:drop, :A), Action(:drop, :B)])])
    @test fix_combinations(d) == Set([Set([Action(:compat, :C)])])
    # ... and every repair it left out is bigger than the ones it offers
    @test all(length(r) > 1 for r in setdiff(repairs, fix_combinations(d)))
    @test d.others === :larger
    report = sprint(show, MIME("text/plain"), d)
    # the menus reach every cheapest repair and only larger ones lie outside,
    # so this is the only *minimal* fix -- not the only fix there is
    @test occursin("The only minimal fix: relax your compat on C", report)
    @test occursin("Larger solutions also exist.", report)
end

@testset "diagnosis: a costlier repair is asked about, not counted" begin
    # Two queries whose menus have the same shape — a product of choices
    # reaching every smallest repair — and that differ only in whether anything
    # beyond those repairs exists at all. The oracle for both is every minimal
    # repair drawn from the actions involved, resolved from scratch.

    # relaxing the bound on :C is the cheap way to settle :A with :B, and
    # dropping both of them is the expensive one; :E with :F is settled by
    # dropping either. So the smallest repairs are two, and two more cost three
    # actions each
    prob = Problem([:A, :B, :E, :F]; compat = Dict(:C => Symbol[]))
    d = check_diagnosis(larger_repair, prob)
    pool = [Action(:compat, :C), Action(:drop, :A), Action(:drop, :B),
            Action(:drop, :E), Action(:drop, :F)]
    repairs = minimal_repairs(pkg_info(larger_repair, prob), prob, pool)
    @test repairs == Set(Set(r) for r in (
        [Action(:compat, :C), Action(:drop, :E)],
        [Action(:compat, :C), Action(:drop, :F)],
        [Action(:drop, :A), Action(:drop, :B), Action(:drop, :E)],
        [Action(:drop, :A), Action(:drop, :B), Action(:drop, :F)]))
    smallest = minimum(length, repairs)
    # the menus are exactly the smallest of them ...
    @test fix_combinations(d) == Set(r for r in repairs if length(r) == smallest)
    # ... and what they leave out costs more, which is what the report says
    @test all(length(r) > smallest
              for r in setdiff(repairs, fix_combinations(d)))
    @test d.others === :larger
    @test occursin("Larger solutions also exist.",
        sprint(show, MIME("text/plain"), d))

    # the same shape of menu over a query where every repair is a smallest one:
    # two conflicts with nothing to do but drop a requirement from each
    prob = Problem([:A, :B, :E, :F])
    d = check_diagnosis(two_conflicts, prob)
    pool = [Action(:drop, p) for p in (:A, :B, :E, :F)]
    repairs = minimal_repairs(pkg_info(two_conflicts, prob), prob, pool)
    @test allequal(length(r) for r in repairs)
    @test fix_combinations(d) == repairs
    @test d.others === :none
    @test !occursin("Other", sprint(show, MIME("text/plain"), d))
end

@testset "diagnosis: one reason, every requirement it breaks" begin
    # :A and :B are each unsatisfiable on their own, for the same reason, and
    # one action rescues both. Either of them alone is enough to make that
    # reason a conflict, so a story that stopped at the smallest one would name
    # one of them and leave the other out of the report altogether — the user
    # would relax the bound on :C and never learn the other had been broken
    prob = Problem([:A, :B];
        compat = Dict(:A => [:a2], :B => [:b2], :C => [:c1]))
    d = check_diagnosis(shared_bound, prob)
    c = only(d.conflicts)
    # both requirements, each with the bound that pins it to the version that
    # needs :C, and the bound on :C they share. Each requirement's own story is
    # told in one piece: what the query left of it, then what that needs
    @test c.reqs == [:A, :B]
    @test [fix.actions for fix in c.fixes] == [[Action(:compat, :C)]]
    @test only(c.fixes).solution == Dict(:A => :a2, :B => :b2, :C => :c2)
    # ... and neither of them can be satisfied on its own, which is what makes
    # leaving one out a lie rather than a shortening
    for p in (:A, :B)
        @test resolve(shared_bound, Problem([p];
            compat = Dict(p => [Symbol("$(lowercase(string(p)))2")],
                          :C => [:c1])); diagnose = false) === nothing
    end
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("you require A, and your compat leaves a2", report)
    @test occursin("you require B, and your compat leaves b2", report)
    @test occursin("A and B cannot both be satisfied", report)
    # both of them fail, and naming one would leave the other unexplained: each
    # gets the line saying what it forces, and `:C` is where they meet
    @test occursin("A a2 requires C c2", report)
    @test occursin("B b2 requires C c2", report)
    @test occursin("your compat leaves C c1", report)
    # dropping both requirements repairs it too, and gives up more
    @test d.others === :larger
    @test occursin("The only minimal fix: relax your compat on C", report)
    # the witness names both of them too
    @test occursin("→ allows: A a2, B b2, C c2", report)
end

@testset "diagnosis: a story ends where its own facts are" begin
    # :E is squeezed by :Z one hop away and by :T two hops away, and the
    # smallest explanation of it names :T. The walk has to end there too: a
    # story that stopped at :Z would state a dependency on a package it then
    # says nothing about, and :T's bound would arrive with nothing in the chain
    # to connect it to — this conflict wearing another one's fact
    prob = Problem([:B, :E];
        compat = Dict(:E => [:e1], :T => [:t1], :U => [:u1], :Z => [:z1]))
    d = check_diagnosis(two_squeezes, prob)
    @test length(d.conflicts) == 2
    c = only(x for x in d.conflicts if x.reqs == [:E])
    # :E is squeezed twice over, so there are two proofs of it and either will
    # do. What is asked of whichever comes back is that it be a proof: that
    # every package it says the query emptied is one its own relations reach,
    # so the chain accounts for the facts it states rather than wearing another
    # conflict's. `check_diagnosis` has already had the harder half -- that
    # every fact is true of the universe, and that the chain is unsatisfiable
    # and minimally so
    reached = Set{Symbol}(c.reqs)
    for l in c.lines, q in packages(l.clause)
        length(packages(l.clause)) > 1 && push!(reached, q)
    end
    @test all(l -> length(packages(l.clause)) > 1 ||
                   only(packages(l.clause)) in reached, c.lines)
    # Either squeeze proves it, and both may be told: they are independent
    # reasons :E cannot be had, and a proof that named one would be a claim
    # that fixing it is enough. What is not allowed is half of one -- a route
    # stated with nothing to close it, which the reachability check above is
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("T", report) || occursin("Z", report)
    @test [fix.actions for fix in c.fixes] ==
        [[Action(:compat, :E)], [Action(:drop, :E)]]
end

@testset "diagnosis: cheapest repairs that are not a product" begin
    # The disagreements form a path :C–:B–:A–:D, so the cheapest repairs are
    # {A,B}, {A,C} and {B,D} — three of them, which is not a product of menu
    # sizes. The report offers the largest rectangle in that family and says
    # there is more of it
    prob = Problem([:A, :B, :C, :D])
    d = check_diagnosis(conflict_path, prob)
    info = pkg_info(conflict_path, prob)
    pool = [Action(:drop, p) for p in (:A, :B, :C, :D)]
    repairs = minimal_repairs(info, prob, pool)
    @test repairs == Set(Set([Action(:drop, x), Action(:drop, y)])
                         for (x, y) in ((:A, :B), (:A, :C), (:B, :D)))
    # what the report offers is a proper part of that, all of it cheapest
    offered = fix_combinations(d)
    @test offered ⊊ repairs
    @test length(offered) == 2
    @test d.others === :some
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("Other solutions also exist.", report)
    @test !occursin("larger", report)
    # A menu of one may only claim as much as `others` knows. Dropping :A is
    # the whole of the first conflict's menu, and it is *not* the only way to
    # settle it -- {drop B, drop D} is a repair the menus never reach, and is
    # no larger. So the report says one fix, not the only one.
    @test occursin("One fix: drop requirement A", report)
    @test !occursin("The only", report)
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


# The report is flat: the query's own facts, then what the registry says about
# them. Nothing on the page is derived from anything else on the page, so a
# line cannot say how its two packages reach each other -- and that is the one
# thing `through` is there to buy back.
@testset "diagnosis: a line names the packages its argument went through" begin
    P, V = String, Int
    VS = Dict(p => [1, 2] for p in ("A", "B", "C", "D", "E"))
    # `p@1 requires q@1`, as the clause it is
    dep(p, q) = Clauses.clause([p => literal(2, [1], true; absent = true),
                                q => literal(2, [1])])
    line(c, through...; pivot = nothing) =
        Line{P}(c, P[through...], false, 1, pivot)
    render(lines; reqs = P[]) = sprint() do io
        Diagnostics.print_conflict(io, Conflict{P,V}(reqs, lines, VS,
            Dict{P,Vector{Vector{Symbol}}}(), Fix{P,V}[]))
    end

    # a line an elimination reached through other packages names them
    @test occursin("A 1 requires E 1 (through B, C and D)",
                   render(Line{P}[line(dep("A", "E"), "B", "C", "D")]))
    # ... and one stated as it stands says nothing about a route
    out = render(Line{P}[line(dep("A", "B"))])
    @test occursin("A 1 requires B 1", out)
    @test !occursin("through", out)

    # where the lines leave one package nothing, saying which saves the reader
    # finding the one name every one of them has in common
    other(p, q) = Clauses.clause([p => literal(2, [1], true; absent = true),
                                  q => literal(2, [2])])
    @test occursin("no version of E is all of these",
                   render(Line{P}[line(dep("A", "E"); pivot = "E"),
                                  line(other("B", "E"); pivot = "E")]))
    # ... and where they do not, there is nothing to say
    @test !occursin("all of these",
                    render(Line{P}[line(dep("A", "E")), line(dep("B", "D"))]))
    # ... including where they do name one package but agree about it: two
    # lines that both leave E at 1 leave it something, and saying otherwise
    # would claim more than the page shows
    @test !occursin("all of these",
                    render(Line{P}[line(dep("A", "E"); pivot = "E"),
                                   line(dep("B", "E"); pivot = "E")]))
    # ... nor where there is only one of them to meet
    @test !occursin("all of these", render(Line{P}[line(dep("A", "E"))]))

    # the query's own facts are said first, whatever order the lines come in,
    # and a requirement is said as one -- not twice, once as a clause
    given = Line{P}(Clauses.clause([
        "B" => literal(2, [1]; absent = true)]), P[], true)
    req = Line{P}(Clauses.clause([
        "A" => literal(2, [1, 2])]), P[], true)
    out = render(Line{P}[line(dep("A", "E")), given, req]; reqs = P["A"])
    @test occursin(r"you require A\n.*B 2 cannot.*\n.*A 1 requires E 1"m, out)
    @test length(collect(eachmatch(r"^  • "m, out))) == 3
end

@testset "diagnosis: the report" begin
    # each conflict carries its own menu, and the menus do not multiply: two
    # menus of two, not one of four. There is no closing sentence — every
    # combination of them is a repair and there are no others. The versions
    # shown are of every package the conflict names, :C and :G included: they
    # are named because the story is about them
    d = resolve(two_conflicts, Problem([:A, :B, :E, :F]))
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 2 conflicts, each of which must be fixed:

        Conflict 1: A and B cannot both be satisfied.
          • you require A
          • you require B
          • no version of C is all of these:
              — A requires C v1
              — B requires C v2
          Fix it by any one of:
            1. drop requirement A
               → allows: B v1, C v2
            2. drop requirement B
               → allows: A v1, C v1

        Conflict 2: E and F cannot both be satisfied.
          • you require E
          • you require F
          • no version of G is all of these:
              — E requires G v1
              — F requires G v2
          Fix it by any one of:
            1. drop requirement E
               → allows: F v1, G v2
            2. drop requirement F
               → allows: E v1, G v1
        """
    # the one-line summary counts the ways of repairing the whole query
    @test sprint(show, d) == "Diagnosis: 2 conflicts, 4 fixes"

    # Every kind that excludes a version is named. Which kind took which
    # version is not: there is no surviving range for the reader to place them
    # against, the fix menu names each kind on its own anyway, and the sentence
    # that would say it reads as the opposite of what it means.
    #
    # And every statement is a line of its own: the availability is on the
    # page once, beside what it contradicts, which is what a premise appearing
    # in a report is for
    d = resolve(needs_dep,
        Problem([:A]; compat = Dict(:B => [:w1]), pin = Dict(:B => :w2)))
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict:

        Conflict 1: A cannot be satisfied.
          • you require A
          • your compat and your pin leaves no version of B
          • A requires B
          Fix it by any one of:
            1. relax your compat on B
               → allows: A v1, B w2
            2. drop requirement A
        """
    @test sprint(show, d) == "Diagnosis: 1 conflict, 2 fixes"
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
                fixes += sum(c -> length(c.fixes), diag.conflicts)
                relaxations += count(a.kind !== :drop
                    for c in diag.conflicts for fix in c.fixes
                    for a in fix.actions)
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
    @test any(l -> packages(l.clause) == ["PrettyTables"], c.lines)
    # relaxing either bound is a fix, and so is not requiring DataFrames
    sets = Set(Set(fix.actions) for fix in c.fixes)
    @test Set([Action(:compat, "PrettyTables")]) ∈ sets
    @test Set([Action(:drop, "DataFrames")]) ∈ sets
    for fix in c.fixes
        fix.actions == [Action(:compat, "PrettyTables")] || continue
        @test fix.solution["DataFrames"] ∈ VersionSpec("1")
        @test haskey(fix.solution, "PrettyTables")
    end
    # the middle of the story: that DataFrames is what needs PrettyTables, and
    # which versions of it do. No bound is stated, since the query has left
    # PrettyTables with nothing whatever DataFrames would have taken
    report = sprint(show, MIME("text/plain"), d)
    # the query's own compat is what took every version away, so it is named:
    # "no version of PrettyTables is available" is the other thing that can
    # empty a package, and it is not this
    @test occursin("your compat leaves no version of PrettyTables", report)
    @test occursin("relax your compat on PrettyTables", report)
    @test occursin("requires PrettyTables", report)

    # ... and a query that leaves PrettyTables something DataFrames will not
    # take is where the bound itself is the story
    prob = Problem(["DataFrames"];
        compat = Dict("DataFrames" => VersionSpec("1.4 - 1.7"),
                      "PrettyTables" => VersionSpec("1")))
    d = check_diagnosis(rp, prob)
    @test d isa Diagnosis{String,VersionNumber}
    c2 = only(d.conflicts)
    @test "DataFrames" in keys(c2.versions) && "PrettyTables" in keys(c2.versions)
    # what it will take is PrettyTables 2, which is what the bound rules out:
    # the statement about the two of them leaves only 2.x of PrettyTables
    report2 = sprint(show, MIME("text/plain"), d)
    @test occursin("requires PrettyTables ", report2)
    for cl in (l.clause for l in c2.lines)
        m = cl["PrettyTables"]
        (m === nothing || length(cl.lits) < 2) && continue
        for (i, v) in enumerate(c2.versions["PrettyTables"])
            m[i] && @test v ∈ VersionSpec("2")
        end
    end

    # a bound that differs across the depending package's versions, on real
    # data: Plots has wanted a different RecipesBase over the years, and a
    # query holding it at 1 and RecipesBase at 0.4 is stopped by all of them
    prob = Problem(["Plots"];
        compat = Dict("Plots" => VersionSpec("1"),
                      "RecipesBase" => VersionSpec("0.4")))
    d = check_diagnosis(rp, prob)
    @test d isa Diagnosis{String,VersionNumber}
    c = only(d.conflicts)
    # The versions of Plots disagree about which RecipesBase they want, but
    # RecipesBase is the only package this proof links Plots to, so nothing it
    # states can tell them apart: what is on the page is about those two and
    # nothing else
    @test Set(p for l in c.lines for p in packages(l.clause)) ==
          Set(["Plots", "RecipesBase"])
end

# The pivot theorem does not promise two sides. Three sets can meet pairwise
# and have nothing in all of them -- and on a line that takes a *disconnected*
# one, since three intervals meeting pairwise share a point. So a three-sided
# meet needs a package one of whose bounds has a hole in it, which is why the
# registry has so few and why one is written out here rather than looked for.
@testset "diagnosis: three sides meeting at one package" begin
    # :P has three versions; each requirement leaves a different pair of them,
    # and :C's is the disconnected one
    data = Dict(
        :A => PkgData([:a1], Dict(:a1 => [:P]), Dict(:a1 => Dict(:P => [:p1, :p2]))),
        :B => PkgData([:b1], Dict(:b1 => [:P]), Dict(:b1 => Dict(:P => [:p2, :p3]))),
        :C => PkgData([:c1], Dict(:c1 => [:P]), Dict(:c1 => Dict(:P => [:p1, :p3]))),
        :P => PkgData([:p3, :p2, :p1], DEPS_NONE, COMP_NONE),
    )
    # every pair of them is fine ...
    for (x, y) in ((:A, :B), (:A, :C), (:B, :C))
        @test resolve(data, Problem([x, y]); diagnose = false) !== nothing
    end
    # ... and the three together are not
    d = check_diagnosis(data, Problem([:A, :B, :C]))
    c = only(d.conflicts)
    @test Set(c.reqs) == Set([:A, :B, :C])
    report = sprint(show, MIME("text/plain"), d)
    # each side says what it demands of the package they meet at, and none of
    # them is left to the reader to infer
    @test occursin("A requires P ≤p2", report)
    @test occursin("B requires P ≥p2", report)
    @test occursin("C requires P p1, p3", report)
    @test occursin("no version of P is all of these", report)
    # ... and every requirement has a demand on it: none of the three is told
    # only by what it rules out
    said = Line{Symbol}[l for l in c.lines if !l.given]
    @test length(said) == 3
    for r in (:A, :B, :C)
        @test any(l -> l.clause[r] !== nothing, said)
    end
end
