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
using Resolver.Diagnostics: Diagnostics, Conflict, Fix, Action, Line, Upstream,
    clause_versions, clauses_satisfiable, clause_of, project, action_phrase,
    report_problems
using Resolver.Clauses: Clauses, Clause, packages, isbottom, clause_phrase,
    literal, resolve_on, subsumes
using Resolver.UnsatCores: sat_mcses

@isdefined(ProofCheck) || include(joinpath(@__DIR__, "proof_check.jl"))
using .ProofCheck
using .ProofCheck: claimed_lines

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

# the repairs one conflict offers: its menu's entries, one by one
offers(c::Conflict) = Diagnostics.selections(c)

# ... and the repairs one block of them offers: one entry from each of its
# conflicts' menus in every combination, and every alternative to it besides
block_offers(d::Diagnosis, g::Vector{Int}) = Diagnostics.block_selections(d, g)

# how the page settles a conflict while it shows what a fix elsewhere gets you:
# the first entry of its menu
default(c::Conflict{P,V}) where {P,V} = unique(first(c.fixes).actions)

# every way of choosing one repair per block, as a set of action sets — the
# fixes of the whole query the report offers
fix_combinations(d::Diagnosis) =
    isempty(d.conflicts) ? Set{Set{Action}}() :
    Set(Set(s) for s in Diagnostics.selections(d))

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
        # because none of them does: each is entailed by the registry rather
        # than by its neighbours. The page reads them as a chain, which is an
        # order and a direction chosen when they are printed, never a claim
        # that one was derived from the last. What used to be checked
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

        # ... and it accounts for what the conflict offers. A fix names
        # packages to act on, and a report that offered one it never spoke of
        # would be talking past itself: the reader is told to change something
        # the argument never mentioned. The conflicts are what the proofs
        # answer for; an alternative repartitions repairs and not reasons, and
        # owes the reader menus and witnesses only.
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
            elsewhere = Action{P}[a for j in eachindex(d.conflicts) if j != i
                                    for a in default(d.conflicts[j])]
            for fix in c.fixes
                @test !isempty(fix.actions)
                @test allunique(fix.actions)
                # a fix shows what it gets you when everything else on the
                # page is settled the first way it offers, so that is what is
                # resolved here — the diagnosis answered the relaxation on the
                # universe filtered for the query, this answers it on a
                # universe filtered for the relaxation itself
                actions = unique!(Action{P}[copy(fix.actions); elsewhere])
                @test fix_resolve(info, prob, actions; order, by) ==
                    fix.solution
            end
            # minimal and distinct: within one menu, no entry asks for a
            # strict superset of another's actions, and no two ask for the
            # same things
            sets = [Set(fix.actions) for fix in c.fixes]
            @test allunique(sets)
            for a in sets, b in sets
                @test !(b ⊊ a)
            end
        end

        # ... and an alternative's entries the same way, held out beside the
        # other menus of its own layer and the conflicts it does not replace
        for a in d.alternatives
            outside = Action{P}[x for j in eachindex(d.conflicts)
                                  if j ∉ a.conflicts
                                  for x in default(d.conflicts[j])]
            for (mi, m) in enumerate(a.menus)
                for fix in m
                    @test !isempty(fix.actions)
                    @test allunique(fix.actions)
                    mates = Action{P}[x for (mj, q) in enumerate(a.menus)
                                        if mj != mi for x in first(q).actions]
                    actions = unique!(Action{P}[copy(fix.actions); mates;
                                                outside])
                    @test fix_resolve(info, prob, actions; order, by) ==
                        fix.solution
                end
                sets = [Set(fix.actions) for fix in m]
                @test allunique(sets)
                for x in sets, y in sets
                    @test !(y ⊊ x)
                end
            end
            # every conflict it replaces is one of its block's, and every menu
            # it avoids is one it replaces
            @test a.avoided ⊆ a.conflicts
            @test !isempty(a.conflicts)
        end

        ## blocked entries
        #
        # An entry answers for an action the page makes tempting and no fix
        # anywhere on it takes, so nothing the cover offers on any layer
        # may say that action. Where the entry names a completion, that
        # completion is the rest of a costlier repair — never a repair the page
        # has printed already (Lemma 28), which would make the road it excuses
        # dead weight instead and the sentence the wrong one
        offered = Set{Action{P}}(a for c in d.conflicts for fix in c.fixes
                                 for a in fix.actions)
        union!(offered, Set{Action{P}}(a for x in d.alternatives
                                         for m in x.menus for fix in m
                                         for a in fix.actions))
        groups = Diagnostics.conflict_blocks(d)
        for c in d.conflicts, (tried, unless) in c.blocks
            @test !isempty(tried)
            @test all(a -> a ∉ offered, tried)
            isempty(unless) && continue
            inside(s) = Set(s) ⊆ Set(unless)
            @test !all(g -> any(inside, block_offers(d, g)), groups)
        end

        ## the cover
        #
        # Each block presents its own share of the cheapest repairs exactly
        # and once: no two of its selections ask for the same things, and none
        # of them is inside another
        for g in groups
            sels = [Set(s) for s in block_offers(d, g)]
            @test allunique(sels)
            for a in sels, b in sels
                @test a == b || !(b ⊊ a)
            end
        end

        ## every combination is a fix
        #
        # The claim the whole report rests on: the conflicts are independent,
        # so any one repair from each of them repairs the query. Checked from
        # scratch, exhaustively — the offers are small enough that they can be
        combos = fix_combinations(d)
        if !isempty(combos) && length(combos) ≤ 64
            for actions in combos
                @test fix_resolve(info, prob, collect(actions);
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

# Two pairs of requirements that press on the same two packages, so their
# repairs entangle: :A, :L and :T all want :C, which the query has narrowed,
# and :L and :T disagree about :D besides. The cheapest repairs are then
# {compat C} × {L, T, compat T} together with {A} × {L, T} — five of them,
# and no product, since dropping :A leaves the reason :L and :T have for
# needing incompatible versions of :C standing unless one of them goes too.
const entangled = Dict(
    :A => PkgData([:a1], Dict(:a1 => [:C]), Dict(:a1 => Dict(:C => [:c3]))),
    :L => PkgData([:l1], Dict(:l1 => [:C, :D]),
        Dict(:l1 => Dict(:C => [:c2, :c3], :D => [:d1]))),
    # :T's versions take disjoint windows of :D, so the prepared universe keeps
    # both of them and relaxing the compat on :T is a fix the menu can offer
    :T => PkgData([:t2, :t1], Dict(:t2 => [:C, :D], :t1 => [:C, :D]),
        Dict(:t2 => Dict(:C => [:c1, :c3], :D => [:d2]),
             :t1 => Dict(:C => [:c1, :c3], :D => [:d1]))),
    :C => PkgData([:c3, :c2, :c1], DEPS_NONE, COMP_NONE),
    :D => PkgData([:d2, :d1], DEPS_NONE, COMP_NONE),
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
    # :R will take. The heading has said the requirement, so the chain opens
    # with what it forces and ends at the fact that contradicts it
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict:

        Conflict 1: R
          • R requires P p2
          • your compat leaves P p1
          Fix it by any one of:
            1. relax your compat on P
               → allows: P p2, R r1
            2. drop requirement R
          Upstream fix: a release of R supporting P p1 would fix this; r1, its latest,
            supports only p2.
            → would allow: P p1
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

        Conflict 1: A
          • A requires C ≥c2
          • your compat leaves C c1
          Fix it by any one of:
            1. relax your compat on C
               → allows: A a3, C c3
            2. drop requirement A
          Upstream fix: a release of A supporting C c1 would fix this; a3, its latest,
            supports only c3.
            → would allow: C c1
        """
end

@testset "diagnosis: a package's availability is a premise where it speaks" begin
    # :C is the package the query emptied and the package the story ends at,
    # and it also speaks: what is left of it needs :A back. The chain says what
    # the query left of it where the package arrives — after the line that asks
    # for it — so the fact stands beside what it contradicts rather than ahead
    # of anything that has introduced :C
    d = check_diagnosis(late_speaker, Problem([:A]; compat = Dict(:C => [:c1])))
    c = only(d.conflicts)
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict:

        Conflict 1: A
          • A requires C c2
          • your compat leaves C c1
          Fix it by any one of:
            1. relax your compat on C
               → allows: A a1, C c2
            2. drop requirement A
          Upstream fix: a release of A supporting C c1 would fix this; a1, its latest,
            supports only c2.
            → would allow: C c1
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
    # is what that sentence is answering, and it is named in full — after the
    # package, which the line names whether or not the heading does
    d2 = check_diagnosis(narrowed_run,
        Problem([:A]; compat = Dict(:A => [:a2, :a1], :B => [:b1])))
    report2 = sprint(show, MIME("text/plain"), d2)
    @test occursin("your compat leaves A ≤a2", report2)
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
    # the chain says :P's bound on :W, then what the query left of :W, and
    # closes with :S's dependency read the other way round -- the same
    # statement, since a clause has no direction, and the way round the chain
    # arrives at it
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("P constrains W w1", report)
    @test occursin("your compat leaves W w2", report)
    @test occursin("W absent leaves no version of S", report)
    # the only dependency stated is the one the registry has: :P's bound on
    # :W permits :W's absence, so nothing on the page says :P brings it in
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
    # the page's own line makes dropping :A tempting and no fix takes it, so
    # the page answers for it: dropping :A alone settles nothing, but it does
    # lie in the costlier minimal repair {:A, :B}, so the verdict names its
    # price rather than dismissing it
    @test occursin("Blocked fixes:", report)
    @test occursin("dropping requirement A would not help unless you also " *
                   "dropped requirement B.", replace(report, r"\n\s+" => " "))
    # that entry is itself the costlier fix, named, so the abstract footer has
    # nothing left to add
    @test !occursin("Costlier fixes also exist.", report)
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
    # what costs more is named where it is tempting: dropping :A lies in a
    # costlier minimal repair, so the verdict states its price, and being a
    # named costlier fix it leaves the abstract footer nothing to say
    report = replace(sprint(show, MIME("text/plain"), d), r"\n\s+" => " ")
    @test occursin("dropping requirement A would not help unless you also " *
                   "dropped requirement B and dropped requirement F.", report)
    @test !occursin("Costlier fixes also exist.", report)

    # the same shape of menu over a query where every repair is a smallest one:
    # two conflicts with nothing to do but drop a requirement from each
    prob = Problem([:A, :B, :E, :F])
    d = check_diagnosis(two_conflicts, prob)
    pool = [Action(:drop, p) for p in (:A, :B, :E, :F)]
    repairs = minimal_repairs(pkg_info(two_conflicts, prob), prob, pool)
    @test allequal(length(r) for r in repairs)
    @test fix_combinations(d) == repairs
    @test d.others === :none
    # every repair is on the menus and nothing costs more, so there is no
    # footer: nothing is left for the page to confess
    report = sprint(show, MIME("text/plain"), d)
    @test !occursin("also exist", report)
    @test !occursin("more minimal fixes", report)
end

@testset "diagnosis: one reason, and what the roads round it cost" begin
    # :A and :B are each unsatisfiable on their own, for the same reason, and
    # one action rescues both. The conflict is about the reason it owns —
    # :A's — and :B is not left out of the report: every road round that
    # reason runs through :B, which is what the blocked entries say
    prob = Problem([:A, :B];
        compat = Dict(:A => [:a2], :B => [:b2], :C => [:c1]))
    d = check_diagnosis(shared_bound, prob)
    c = only(d.conflicts)
    # the conflict answers for the requirement its own reason uses, with the
    # bound that pins it to the version that needs :C and the bound on :C it
    # collides with. That story is told in one piece
    @test c.reqs == [:A]
    @test [fix.actions for fix in c.fixes] == [[Action(:compat, :C)]]
    @test only(c.fixes).solution == Dict(:A => :a2, :B => :b2, :C => :c2)
    # ... and neither requirement can be satisfied on its own, which is why no
    # action on :A alone is on the menu
    for p in (:A, :B)
        @test resolve(shared_bound, Problem([p];
            compat = Dict(p => [Symbol("$(lowercase(string(p)))2")],
                          :C => [:c1])); diagnose = false) === nothing
    end
    report = sprint(show, MIME("text/plain"), d)
    wrapped = replace(report, r"\n\s+" => " ")
    # the heading names the requirement and claims nothing about it: the
    # conflict reads under an implicit "given the rest of the requirements"
    @test occursin("Conflict 1: A\n", report)
    # the page tells that reason in full: what the query left of :A, what that
    # forces, and the bound on :C it collides with
    @test occursin("your compat leaves A a2", report)
    @test occursin("A a2 requires C c2", report)
    @test occursin("your compat leaves C c1", report)
    # the two actions on :A the page makes tempting and no fix takes, each
    # with the solver's verdict: each lies in a costlier minimal repair — the
    # bound on :A with :B given up, the requirement on :A with :B's bound
    # relaxed — which is where :B is accounted for. A tried bundle is a vector
    # per tempting fact, so each entry's first component is nested
    @test c.blocks == [([[Action(:compat, :A)]], [Action(:drop, :B)]),
                       ([[Action(:drop, :A)]], [Action(:compat, :B)])]
    @test occursin("relaxing your compat on A would not help unless you " *
                   "also dropped requirement B.", wrapped)
    @test occursin("dropping requirement A would not help unless you also " *
                   "relaxed your compat on B.", wrapped)
    # :B's own reason proves the same conflict a second way and prints
    # nothing: one conflict, one story, and no line of it names :B
    @test !occursin("your compat leaves B b2", report)
    @test !occursin("B b2 requires C c2", report)
    @test !any(l -> :B in packages(l.clause), c.lines)
    # dropping both requirements repairs it too, and gives up more -- which
    # the unless entry above has already named concretely, so the abstract
    # announcement is left off while the fact behind it stands
    @test d.others === :larger
    @test !occursin("Costlier fixes also exist.", report)
    @test occursin("The only minimal fix: relax your compat on C", report)
    # the witness is shown of the packages the page speaks of
    @test occursin("→ allows: A a2, C c2", report)
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
    # {A,B}, {A,C} and {B,D} — three of them, which is no product of menu
    # sizes. The largest rectangle in the family leads and its menus are the
    # conflicts, and what it does not reach follows them as the other way it is
    prob = Problem([:A, :B, :C, :D])
    d = check_diagnosis(conflict_path, prob)
    info = pkg_info(conflict_path, prob)
    pool = [Action(:drop, p) for p in (:A, :B, :C, :D)]
    repairs = minimal_repairs(info, prob, pool)
    @test repairs == Set(Set([Action(:drop, x), Action(:drop, y)])
                         for (x, y) in ((:A, :B), (:A, :C), (:B, :D)))
    # the family does not factor, so its leading layer's two menus are the two
    # conflicts and what that product does not reach is the alternative to
    # them both: nothing offered that is not a cheapest repair, and no cheapest
    # repair left out
    @test length(d.conflicts) == 2
    @test [[[a.pkg for a in f.actions] for f in c.fixes] for c in d.conflicts] ==
        [[[:A]], [[:B], [:C]]]
    a = only(d.alternatives)
    @test a.conflicts == [1, 2]
    # two menus of one entry each, not one compound entry: the split is
    # already exact, so nothing couples them, and they print as one line
    @test [[[x.pkg for x in f.actions] for f in m] for m in a.menus] ==
        [[[:B]], [[:D]]]
    @test fix_combinations(d) == repairs
    # every selection of the alternative misses the whole of the first menu,
    # and that is what the label says it is doing without
    @test a.avoided == [1]
    # so nothing is outside the cover, and the enumeration was not cut short
    @test d.others === :none
    report = sprint(show, MIME("text/plain"), d)
    # the alternative denies that the solutions are one fix from each menu, so
    # the headline does not claim it
    @test startswith(report, "Unsatisfiable — 2 conflicts, pick a fix for each:")
    # each conflict argues the one reason its own menu owns, as a chain of its
    # own: pooled into one they would read as a single argument that is none
    @test all(length(unique(l.proof for l in c.lines)) == 1 for c in d.conflicts)
    @test !occursin("  and also:\n", report)
    wrapped = replace(report, r"\n\s+" => " ")
    @test occursin("One fix: drop requirement A", wrapped)
    @test occursin("1. drop requirement B", wrapped)
    @test occursin("2. drop requirement C", wrapped)
    @test occursin("Or, to fix without dropping requirement A: drop " *
                   "requirement B and drop requirement D", wrapped)
    # ... and the alternative states fixes only: reasons do not layer
    tail = split(report, "Or, to fix without")[2]
    @test !occursin("requires", tail)
    @test !occursin("Blocked fixes", tail)
    # nothing is left over, so the page has no residue to announce
    @test !occursin("If none of the fixes above suits", report)
    @test !occursin("also exist", report)
    @test !occursin("more minimal fixes than are shown", report)
end

@testset "diagnosis: a cover completes an entangled family" begin
    # :A, :L and :T all want :C, which the query has narrowed, and :L and :T
    # disagree about :D besides. The cheapest repairs are five and no product
    # of anything: a 2×2 rectangle leads — {compat C, drop A} against
    # {drop L, drop T}, one conflict apiece — and the fifth repair, which no
    # rectangle of two reaches beside them, is the alternative to them both
    prob = Problem([:A, :L, :T]; compat = Dict(:C => [:c1, :c2], :T => [:t2]))
    d = check_diagnosis(entangled, prob)
    info = pkg_info(entangled, prob)
    pool = [Action(:compat, :C), Action(:compat, :T), Action(:drop, :A),
            Action(:drop, :L), Action(:drop, :T)]
    repairs = minimal_repairs(info, prob, pool)
    @test length(repairs) == 5
    @test all(length(r) == 2 for r in repairs)
    @test length(d.conflicts) == 2
    a = only(d.alternatives)
    # the page reaches every cheapest repair there is, and offers nothing else
    @test fix_combinations(d) == repairs
    @test [[Set(f.actions) for f in c.fixes] for c in d.conflicts] ==
        [[Set([Action(:compat, :C)]), Set([Action(:drop, :A)])],
         [Set([Action(:drop, :L)]), Set([Action(:drop, :T)])]]
    @test a.conflicts == [1, 2]
    # two menus of one entry each: the split is exact, so nothing couples
    # them, and the page says them as one line
    @test [[Set(f.actions) for f in m] for m in a.menus] ==
        [[Set([Action(:compat, :C)])], [Set([Action(:compat, :T)])]]
    # each entry carries a witness of its own, taken with its layer-mates at
    # their first entries, which resolves
    for (j, m) in enumerate(a.menus), fix in m
        @test !isempty(fix.solution)
        mates = [first(a.menus[k]).actions for k in eachindex(a.menus) if k != j]
        @test fix_resolve(info, prob, vcat(fix.actions, mates...)) == fix.solution
    end
    report = sprint(show, MIME("text/plain"), d)
    @test startswith(report, "Unsatisfiable — 2 conflicts, pick a fix for each:")
    wrapped = replace(report, r"\n\s+" => " ")
    @test occursin("1. relax your compat on C", wrapped)
    @test occursin("2. drop requirement A", wrapped)
    @test occursin("1. drop requirement L", wrapped)
    @test occursin("2. drop requirement T", wrapped)
    # It relaxes the compat on C, which the first menu also offers, so there is
    # no single thing it does without there — what it declines whole is the
    # second menu, and the label names that conflict rather than a fix
    @test a.avoided == [2]
    @test occursin("Or, to fix without any of the fixes for Conflict 2: " *
                   "relax your compat on C and relax your compat on T → " *
                   "allows: A a1, C c3, L l1, T t1", wrapped)
    # and no proof under the alternative: reasons do not layer, so the
    # conflicts above have already explained every one there is and what
    # follows states fixes only
    tail = split(report, "Or, to fix without any of the fixes for")[2]
    @test !occursin("requires", tail)
    @test !occursin("your compat leaves", tail)
    @test !occursin("Blocked fixes", tail)
    # nothing is outside the cover and nothing costlier exists, so the page
    # has no gap to confess and prints no footer
    @test d.others === :none
    @test !occursin("also exist", report)
    @test !occursin("more minimal fixes than are shown", report)
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


# The report is a chain: a root fact, the statements that carry it, and the
# fact it ends against. A line still cannot say how its two packages reach each
# other -- the elimination between them is not on the page -- and that is the
# one thing `through` is there to buy back.
@testset "diagnosis: a line names the packages its argument went through" begin
    P, V = String, Int
    VS = Dict(p => [1, 2] for p in ("A", "B", "C", "D", "E"))
    # `p@1 requires q@1`, as the clause it is
    dep(p, q) = Clauses.clause([p => literal(2, [1], true; absent = true),
                                q => literal(2, [1])])
    line(c, through...; pivot = nothing) =
        Line{P}(c, P[through...], false, 1, pivot)
    render(lines; reqs = P[], vs = VS) = sprint() do io
        Diagnostics.print_conflict(io, Conflict{P,V}(reqs, lines, vs,
            Dict{P,Vector{Vector{Symbol}}}(), Fix{P,V}[]))
    end

    # a line an elimination reached through other packages names them
    @test occursin("A 1 requires E 1 (through B, C and D)",
                   render(Line{P}[line(dep("A", "E"), "B", "C", "D")]))
    # ... and one stated as it stands says nothing about a route
    out = render(Line{P}[line(dep("A", "B"))])
    @test occursin("A 1 requires B 1", out)
    @test !occursin("through", out)

    # where three lines leave one package nothing, saying which saves the
    # reader finding the one name every one of them has in common
    VS3 = Dict(p => [1, 2, 3] for p in ("A", "B", "C", "E"))
    at(p, vs) = Clauses.clause([p => literal(3, [1], true; absent = true),
                                "E" => literal(3, vs)])
    @test occursin("incompatible constraints on E:",
                   render(Line{P}[line(at("A", [1, 2]); pivot = "E"),
                                  line(at("B", [2, 3]); pivot = "E"),
                                  line(at("C", [1, 3]); pivot = "E")];
                          vs = VS3))
    # ... and two sides of one package are a chain rather than a meet: a clause
    # has no direction, so the second is said to the package it rests on and
    # read against what the first left. (This fixture was the two-sided meet
    # display; two sides linearize always, so the display is gone from it.)
    other(p, q) = Clauses.clause([p => literal(2, [1], true; absent = true),
                                  q => literal(2, [2])])
    two = render(Line{P}[line(dep("A", "E"); pivot = "E"),
                         line(other("B", "E"); pivot = "E")])
    @test !occursin("incompatible constraints", two)
    @test occursin("A 1 requires E 1", two)
    @test occursin("E 1 constrains B 2", two)
    # ... and where the lines do not leave the package nothing, there is
    # nothing to say
    @test !occursin("incompatible constraints",
                    render(Line{P}[line(dep("A", "E")), line(dep("B", "D"))]))
    # ... including where they do name one package but agree about it: lines
    # that all leave E at 1 leave it something, and saying otherwise would
    # claim more than the page shows
    @test !occursin("incompatible constraints",
                    render(Line{P}[line(at("A", [1]); pivot = "E"),
                                   line(at("B", [1]); pivot = "E"),
                                   line(at("C", [1]); pivot = "E")];
                           vs = VS3))
    # ... nor where there is only one of them to meet
    @test !occursin("incompatible constraints",
                    render(Line{P}[line(dep("A", "E"))]))

    # a limit on a package no statement reaches ends the page rather than
    # opening it -- the chain is what the reader is following -- and the
    # requirement is not said at all: the heading names it, so a line for it
    # would be the page saying one thing twice
    given = Line{P}(Clauses.clause([
        "B" => literal(2, [1]; absent = true)]), P[], true, 1)
    req = Line{P}(Clauses.clause([
        "A" => literal(2, [1, 2])]), P[], true, 1)
    out = render(Line{P}[line(dep("A", "E")), given, req]; reqs = P["A"])
    @test occursin(r"A 1 requires E 1.*\n.*B 2 cannot"m, out)
    @test !occursin("you require", out)
    @test length(collect(eachmatch(r"^  • "m, out))) == 2
end

# Licensed coarsening: a parallel family joins by resolution on one of its
# packages exactly when the claim still contradicts afterwards. A staircase of
# thresholds the proof never needs collapses to one line; a boundary the
# contradiction stands on refuses to.
@testset "diagnosis: joins are licensed by the claim, not by a local rule" begin
    # A's three versions pick disjoint windows of C — so the prepared universe
    # keeps them all — and B needs C at v4. The clauses under test are built
    # over the instance's version lists; `clauses_satisfiable` is pure logic
    # over those domains, so the registry behind them only has to keep the
    # vocabulary alive.
    data = Dict(
        :A => PkgData([:v3, :v2, :v1],
            Dict(v => [:C] for v in (:v1, :v2, :v3)),
            Dict(:v1 => Dict(:C => [:v1]),
                 :v2 => Dict(:C => [:v2]),
                 :v3 => Dict(:C => [:v3]))),
        :B => PkgData([:v1], Dict(:v1 => [:C]), Dict(:v1 => Dict(:C => [:v4]))),
        :C => PkgData([:v4, :v3, :v2, :v1], Dict{Symbol,Vector{Symbol}}(),
                      Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()))
    sat, univ, info = failed_instance(data, Problem([:A, :B]))
    vers(p) = clause_versions(sat, p)
    order(p) = Clauses.version_order(vers(p))
    ix(p, v) = findfirst(==(v), vers(p))
    imp(p, R, q, S) = Clauses.clause([
        p => literal(length(vers(p)), Int[ix(p, v) for v in R], true; absent = true),
        q => literal(length(vers(q)), Int[ix(q, v) for v in S])])
    stair = [imp(:A, [:v1], :C, [:v1]),
             imp(:A, [:v2], :C, [:v1, :v2]),
             imp(:A, [:v3], :C, [:v1, :v2, :v3])]
    held = [Clauses.clause([:A => literal(3, 1:3)]),
            Clauses.clause([:B => literal(1, 1:1)]),
            imp(:B, [:v1], :C, [:v4])]
    # the staircase's thresholds never matter: every A caps C below v4, so the
    # whole family joins to one clause and the claim still contradicts
    out = Diagnostics.coarsen_core(sat, Vector{Clause{Symbol}}(stair),
                                   Vector{Clause{Symbol}}(held))
    @test length(out) == 1
    @test !clauses_satisfiable(sat, [held; out])
    # ... but a boundary the contradiction stands on refuses to join: here B
    # tolerates C v3, so "A v3 caps C at v3" is the one line that closes, and
    # joining it away would leave the claim satisfiable
    held2 = [Clauses.clause([:A => literal(3, [ix(:A, :v3)])]),
             Clauses.clause([:B => literal(1, 1:1)]),
             imp(:B, [:v1], :C, [:v3, :v4]),
             Clauses.clause([:C => literal(4, [ix(:C, :v3)], true; absent = true)])]
    out2 = Diagnostics.coarsen_core(sat, Vector{Clause{Symbol}}(stair),
                                    Vector{Clause{Symbol}}(held2))
    @test !clauses_satisfiable(sat, [held2; out2])
end

# A conflict's further reasons print after the menu, as blocked fixes: each
# holds with its entry's actions withdrawn, so what it argues is that those
# actions settle nothing. One sentence per entry — the actions in the trying
# and their verdict — with the proof behind it, in the lines, unprinted.
# A conflict whose share of the repairs is not one choose-one menu prints as
# the layers it is covered by: the leading product as bullets to settle, and
# what it does not reach after it. Built by hand, since what is under test is
# the layout and not the analysis that found the family.
@testset "diagnosis: an alternative prints after the conflicts" begin
    P, V = String, Int
    VS = Dict("X" => [1, 2])
    lines = Line{P}[Line{P}(Clauses.clause(["X" => literal(2, [1])]), P[], true)]
    fix(as...) = Fix{P,V}(Action{P}[Action(:drop, a) for a in as],
                          Dict("X" => 1))
    conflict(fs...) = Conflict{P,V}(P["X"], lines, VS,
        Dict{P,Vector{Vector{Symbol}}}(), Fix{P,V}[fs...])
    menu(fs...) = Fix{P,V}[fs...]
    alt(cs, av, ms...) = Diagnostics.Alternative{P,V}(
        Int[cs...], Int[av...], Vector{Fix{P,V}}[ms...])
    page(cs, as) = sprint(show, MIME("text/plain"),
                          Diagnosis(Conflict{P,V}[cs...],
                                    Diagnostics.Alternative{P,V}[as...], :none))

    # An alternative to a block of two conflicts, one of which it declines
    # whole: the label says the fix it is doing without, and what it offers
    # instead is its one-entry menus said as one thing to do, then the menu
    # with a choice numbered like a conflict's, one witness under each entry
    two = page((conflict(fix("A")), conflict(fix("B"), fix("C"))),
               (alt([1, 2], [1], menu(fix("D", "E")), menu(fix("F"), fix("G"))),))
    @test occursin("Or, to fix without dropping requirement A:\n" *
                   "  drop requirement D and drop requirement E, and one of:\n" *
                   "    1. drop requirement F\n" *
                   "       → allows: X 1\n" *
                   "    2. drop requirement G\n" *
                   "       → allows: X 1\n", two)
    @test !occursin("•", split(two, "Or, to fix")[2])
    # ... and it prints after the last conflict, not inside either of them
    @test findfirst("Or, to fix", two)[1] > findlast("Conflict 2", two)[1]
    # a layer that is one menu with a choice, and nothing else, is that choice
    pick = page((conflict(fix("A")), conflict(fix("B"))),
                (alt([1, 2], [1], menu(fix("F"), fix("G"))),))
    @test occursin("Or, to fix without dropping requirement A:\n  any one of:\n" *
                   "    1. drop requirement F\n", pick)
    # ... and several menus with a choice are settled each, as bullets
    several = page((conflict(fix("A")), conflict(fix("B"))),
                   (alt([1, 2], [1], menu(fix("D")), menu(fix("F"), fix("G")),
                        menu(fix("H"), fix("I"))),))
    @test occursin("  drop requirement D, and settle each of these:\n" *
                   "  • drop requirement F, or drop requirement G\n", several)
    @test occursin("  • drop requirement H, or drop requirement I\n", several)

    # an alternative that leaves no choice at all is one thing to do, and
    # prints as the line it is
    one = page((conflict(fix("A")), conflict(fix("B"))),
               (alt([1, 2], [1, 2], menu(fix("C", "D"))),))
    @test occursin("Or, to fix without dropping requirement A or dropping " *
                   "requirement B:\n  drop requirement C and drop " *
                   "requirement D\n  → allows: X 1", one)
    # ... and so does one whose several menus each hold a single entry: they
    # are all to be settled, and bulleted they would read as the choice they
    # are not
    both = page((conflict(fix("A")), conflict(fix("B"))),
                (alt([1, 2], [1, 2], menu(fix("C")), menu(fix("D"))),))
    @test occursin("  drop requirement C and drop requirement D\n" *
                   "  → allows: X 1", both)
    @test !occursin("•", split(both, "Or, to fix")[2])

    # where an avoided menu offers a choice there is no single fix to name, so
    # the label points at the conflicts it replaces instead
    many = page((conflict(fix("A"), fix("B")), conflict(fix("C")),
                 conflict(fix("D"), fix("E"))),
                (alt([1, 2, 3], [1, 3], menu(fix("F"))),))
    @test occursin("Or, to fix without any of the fixes for Conflicts 1 and 3:", many)
    one_of = page((conflict(fix("A"), fix("B")), conflict(fix("C"))),
                  (alt([1, 2], [1], menu(fix("F"))),))
    @test occursin("Or, to fix without any of the fixes for Conflict 1:", one_of)

    # an alternative every conflict of its block has a fix inside declines
    # nothing whole, and says only that it is another way
    none = page((conflict(fix("A"), fix("B")), conflict(fix("C"))),
                (alt([1, 2], Int[], menu(fix("A", "C"))),))
    @test occursin("Or, to fix another way:", none)

    # a conflict of a block with no alternative is the whole of what settles
    # that block, and its menu of one may say so; one whose block has an
    # alternative may not
    plain = page((conflict(fix("A")), conflict(fix("B"), fix("C"))), ())
    @test occursin("The only fix: drop requirement A", plain)
    @test occursin("Fix it by any one of:", plain)
    @test !occursin("Or,", plain)
    # ... and the claim is block-local: the first conflict here is settled by
    # its own menu whatever the second block's alternative offers
    split_blocks = page((conflict(fix("A")), conflict(fix("B"))),
                        (alt([2], [2], menu(fix("C", "D"))),))
    @test occursin("The only fix: drop requirement A", split_blocks)
    @test occursin("One fix: drop requirement B", split_blocks)
end

@testset "diagnosis: the headline tells the reader what to do" begin
    P, V = String, Int
    VS = Dict("X" => [1, 2])
    lines = Line{P}[Line{P}(Clauses.clause(["X" => literal(2, [1])]), P[], true)]
    fix(as...) = Fix{P,V}(Action{P}[Action(:drop, a) for a in as],
                          Dict("X" => 1))
    conflict(fs...) = Conflict{P,V}(P["X"], lines, VS,
        Dict{P,Vector{Vector{Symbol}}}(), Fix{P,V}[fs...])
    head(cs, as) = first(split(sprint(show, MIME("text/plain"),
        Diagnosis(Conflict{P,V}[cs...],
                  Diagnostics.Alternative{P,V}[as...], :none)), "\n"))

    # With several conflicts the headline is an instruction: one fix from
    # each menu, in every combination, is a cheapest repair
    @test head((conflict(fix("A")), conflict(fix("B"), fix("C"))), ()) ==
        "Unsatisfiable — 2 conflicts, pick a fix for each:"
    # An alternative is a cheapest repair taking no entry of some menu; the
    # instruction stands, and the alternative says the other way itself,
    # opening with "Or," — so the two are never read as one claim
    alt = Diagnostics.Alternative{P,V}([1, 2], [1],
                                       Vector{Fix{P,V}}[Fix{P,V}[fix("D")]])
    @test head((conflict(fix("A")), conflict(fix("B"), fix("C"))), (alt,)) ==
        "Unsatisfiable — 2 conflicts, pick a fix for each:"
    @test occursin("\nOr, to fix without dropping requirement A:\n",
        sprint(show, MIME("text/plain"),
               Diagnosis(Conflict{P,V}[conflict(fix("A")), conflict(fix("B"), fix("C"))],
                         Diagnostics.Alternative{P,V}[alt], :none)))
    # one conflict has nothing to pick between
    @test head((conflict(fix("A")),), ()) == "Unsatisfiable — 1 conflict:"
end

@testset "diagnosis: blocked fixes print after the menu" begin
    # The section is indexed by action, not by reason: one sentence for each
    # action the page makes tempting and no fix takes, and the verdict in it
    # was decided by the solver long before this printer saw it. No proof
    # prints here -- why a fix is not offered is a second-order question, and
    # the verdict has already answered it.
    P, V = String, Int
    VS = Dict(p => [1, 2] for p in ("A", "B", "E"))
    dep(p, q, v) = Clauses.clause([p => literal(2, [1], true; absent = true),
                                   q => literal(2, [v])])
    lines = Line{P}[Line{P}(dep("A", "B", 1), P[], false, 1, nothing)]
    entries(bs...) = Tuple{Vector{Vector{Action{P}}},Vector{Action{P}}}[bs...]
    page(blocks, index = nothing) = sprint() do io
        Diagnostics.print_conflict(io, Conflict{P,V}(P["A", "B"], lines, VS,
            Dict{P,Vector{Vector{Symbol}}}(), Fix{P,V}[], blocks), index)
    end

    idle = entries(([[Action(:compat, "A")]], Action{P}[]))
    out = page(idle)
    @test occursin("Blocked fixes:", out)
    @test occursin("relaxing your compat on A does not help.", out)
    # the body stays above it: the reader meets the argument and the offer
    # first, and the roads not taken second
    @test first(findfirst("A 1 requires B 1", out)) <
          first(findfirst("Blocked fixes:", out))

    # a completion turns the flat refusal into what the road would cost
    outu = page(entries(([[Action(:compat, "A")]], [Action(:drop, "B")])))
    @test occursin("relaxing your compat on A would not help unless you " *
                   "also dropped requirement B.",
                   replace(outu, r"\n\s+" => " "))

    # ... named in full while it is short enough to act on, and counted once
    # it is a verdict rather than a list: past four further actions the
    # sentence says how many, not which
    long = [Action(:drop, p) for p in ("B", "C", "D", "E", "F")]
    outl = page(entries(([[Action(:compat, "A")]], long)))
    @test occursin("relaxing your compat on A would not help without 5 " *
                   "other changes.", replace(outl, r"\n\s+" => " "))
    @test !occursin("dropped requirement B", outl)
    outf = page(entries(([[Action(:compat, "A")]], long[1:4])))
    @test occursin("would not help unless you also dropped requirement B, " *
                   "dropped requirement C, dropped requirement D and dropped " *
                   "requirement E.", replace(outf, r"\n\s+" => " "))
    outgl = page(entries(([[Action(:compat, "A")], [Action(:drop, "G")]], long)))
    @test occursin("would only help if you do both and 5 other changes.",
                   replace(outgl, r"\n\s+" => " "))

    # two tempting actions that exhibit one repair are one entry: each is
    # insufficient alone and the repair carries them both, so the page says
    # the repair once from both its ends
    outg = page(entries(([[Action(:compat, "A")], [Action(:drop, "B")]],
                         Action{P}[])))
    @test occursin("relaxing your compat on A or dropping requirement B " *
                   "would only help if you do both.",
                   replace(outg, r"\n\s+" => " "))
    # lifting several kinds is one edit of the reader's, so it is one bundle,
    # and the verb agrees with the actions rather than with the entry
    outm = page(entries(([[Action(:compat, "A"), Action(:pin, "A")]],
                         Action{P}[])))
    @test occursin("relaxing your compat on A and unpinning A do not help.",
                   replace(outm, r"\n\s+" => " "))

    # the heading is the requirements the conflict answers for, bare: it says
    # what the conflict is about and claims nothing about them
    @test occursin("Conflict 1: A and B\n", page(idle, 1))

    # a conflict with nothing tempting says nothing about blocked fixes
    @test !occursin("Blocked fixes", page(entries()))
end


# A heading is the requirements its conflict answers for, so two conflicts
# rooted in the same requirement print the same one -- which reads as one
# conflict said twice rather than as two problems sharing a root. Where
# headings collide, each is extended with the package that conflict closes
# against; a heading nothing collides with is untouched, and conflicts whose
# lines say the same thing leave the number to do the whole of the telling.
# Built by hand, since what is under test is the rendering and not the
# analysis that found the conflicts.
@testset "diagnosis: a repeated heading names what tells the two apart" begin
    P, V = String, Int
    VS = Dict(p => [1, 2] for p in ("A", "B", "C", "D", "E"))
    dep(p, q) = Clauses.clause([p => literal(2, [1], true; absent = true),
                                q => literal(2, [1])])
    bound(q) = Clauses.clause([q => literal(2, [2]; absent = true)])
    # `p 1 requires q 1`, and the user's compat on q contradicting it
    conflict(p, q) = Conflict{P,V}(P[p],
        Line{P}[Line{P}(dep(p, q), P[], false, 1, nothing),
                Line{P}(bound(q), P[], true)],
        VS, Dict{P,Vector{Vector{Symbol}}}(q => [[:compat]]), Fix{P,V}[])
    page(cs...) = replace(sprint(show, MIME("text/plain"),
                                 Diagnosis(Conflict{P,V}[cs...], :none)),
                          r"\n\s+" => " ")

    # two conflicts about A closing at different packages: each heading names
    # the package that conflict contradicts at, and neither names the other's
    out = page(conflict("A", "B"), conflict("A", "C"), conflict("D", "E"))
    @test occursin("Conflict 1: A and B • A 1 requires B 1", out)
    @test occursin("Conflict 2: A and C • A 1 requires C 1", out)
    # ... and the conflict nothing collides with keeps its bare heading
    @test occursin("Conflict 3: D • D 1 requires E 1", out)

    # a heading with nothing to collide with is untouched, same packages or no
    @test occursin("Conflict 1: A • A 1 requires B 1",
                   page(conflict("A", "B")))

    # two conflicts whose lines say the same thing have nothing to add: the
    # number is the whole of the difference
    same = page(conflict("A", "B"), conflict("A", "B"))
    @test occursin("Conflict 1: A • A 1 requires B 1", same)
    @test occursin("Conflict 2: A • A 1 requires B 1", same)
    @test !occursin("A and B", same)
end

@testset "diagnosis: the report" begin
    # each conflict carries its own menu, and the menus do not multiply: two
    # menus of two, not one of four. There is no closing sentence — every
    # combination of them is a repair and there are no others. The versions
    # shown are of every package the conflict names, :C and :G included: they
    # are named because the story is about them
    d = resolve(two_conflicts, Problem([:A, :B, :E, :F]))
    # the menus are a genuine product — every combination is a repair and
    # there are no others — so "pick a fix for each" is the whole of it: the
    # solutions really are one-fix-from-each-menu
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 2 conflicts, pick a fix for each:

        Conflict 1: A and B
          • A requires C v1
          • C v1 leaves no version of B
          Fix it by any one of:
            1. drop requirement A
               → allows: B v1, C v2
            2. drop requirement B
               → allows: A v1, C v1

        Conflict 2: E and F
          • E requires G v1
          • G v1 leaves no version of F
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

        Conflict 1: A
          • A requires B
          • your compat and your pin leaves no version of B
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
    @test occursin("incompatible constraints on P:", report)
    # ... and every requirement has a demand on it: none of the three is told
    # only by what it rules out
    said = Line{Symbol}[l for l in c.lines if !l.given]
    @test length(said) == 3
    for r in (:A, :B, :C)
        @test any(l -> l.clause[r] !== nothing, said)
    end
end


# An upstream fix is the one thing on the page the reader cannot do themselves:
# a release of some package, supporting a package the query narrowed, that this
# has resolved and found would settle the conflict. The bar is three conditions
# (Section 7 of the theory page) and every one of them is checked here on data
# small enough to read: the bound has to meet one of the user's *own* facts, the
# query has to admit the releasing package's latest, and the solve has to
# succeed with that latest taken.
#
# `report_problems` is the checker, and it is given the query and the data, so
# the four questions only the registry can answer -- is that the latest, is that
# its bound, does the witness land outside it, does the query admit it -- are
# asked here too.

# :A's only version wants the older :B, and the user wants the newer one: the
# bound is :A's and lifting it is not the user's to do
const upstream_pair = Dict(
    :A => PkgData([:a1], Dict(:a1 => [:B]), Dict(:a1 => Dict(:B => [:b1]))),
    :B => PkgData([:b2, :b1], DEPS_NONE, COMP_NONE),
)

@testset "diagnosis: a release someone else could cut" begin
    prob = Problem([:A, :B]; compat = Dict(:B => [:b2]))
    d = check_diagnosis(upstream_pair, prob)
    c = only(d.conflicts)
    u = only(c.upstream)
    @test (u.pkg, u.latest, u.dep, u.supports) == (:A, :a1, :B, :b2)
    @test u.supported == [:b1]
    # the versions are a resolve's, like every other witness on the page
    @test u.solution == Dict(:A => :a1, :B => :b2)
    # ... and the sentence is one the reader could send as it stands (read
    # here with the wrapping undone, which is the terminal's business)
    report = sprint(show, MIME("text/plain"), d)
    unwrapped = replace(report, r"\n\s+" => " ")
    @test occursin("Upstream fix: a release of A supporting B b2 would fix " *
                   "this; a1, its latest, supports only b1.", unwrapped)
    @test occursin("→ would allow: B b2", report)
    # (V7) everything Section 8 asks of it, the registry's part included
    @test isempty(report_problems(d; prob, data = upstream_pair))
    @test isempty(report_problems(d))
    # ... and the probes are the resolve's to run, so a caller can say no
    d2 = resolve(upstream_pair, prob; upstream = false)
    @test all(c -> isempty(c.upstream), d2.conflicts)
    @test !occursin("Upstream", sprint(show, MIME("text/plain"), d2))
end

@testset "diagnosis: a checker that reads the registry" begin
    # V7 is a check and not a formality: a sentence the data does not bear out
    # is caught, whichever half of it is wrong
    prob = Problem([:A, :B]; compat = Dict(:B => [:b2]))
    d = check_diagnosis(upstream_pair, prob)
    c = only(d.conflicts)
    u = only(c.upstream)
    function retold(v::Upstream{Symbol,Symbol})
        c2 = Conflict{Symbol,Symbol}(c.reqs, c.lines, c.versions, c.excluded,
                                     c.fixes, c.blocks, [v])
        return report_problems(Diagnosis([c2], d.others); prob,
                               data = upstream_pair)
    end
    @test isempty(retold(u))
    # a witness that does not take the release
    @test !isempty(retold(Upstream(:A, :a1, :B, :b2, [:b1],
                                   Dict(:B => :b2))))
    # a version that is not the latest
    @test !isempty(retold(Upstream(:B, :b1, :A, :a1, Symbol[],
                                   Dict(:A => :a1, :B => :b1))))
    # a version of the bounded package the bound admits after all: then the
    # release drops a bound that was not in the way (Lemma 32)
    @test !isempty(retold(Upstream(:A, :a1, :B, :b1, [:b1],
                                   Dict(:A => :a1, :B => :b1))))
    # a package no line of the conflict says the query narrowed
    @test !isempty(retold(Upstream(:B, :b2, :A, :a1, Symbol[],
                                   Dict(:A => :a1, :B => :b2))))
end

@testset "diagnosis: a bound the user's own facts do not meet" begin
    # :A wants the old :C, :B wants the new one, and the query says nothing
    # about :C at all. Two maintainers could each fix this and the page will
    # not judge between them, so it asks neither
    data = Dict(
        :A => PkgData([:a1], Dict(:a1 => [:C]), Dict(:a1 => Dict(:C => [:c1]))),
        :B => PkgData([:b1], Dict(:b1 => [:C]), Dict(:b1 => Dict(:C => [:c2]))),
        :C => PkgData([:c2, :c1], DEPS_NONE, COMP_NONE),
    )
    d = check_diagnosis(data, Problem([:A, :B]))
    @test all(c -> isempty(c.upstream), d.conflicts)
    @test !occursin("Upstream", sprint(show, MIME("text/plain"), d))
    @test isempty(report_problems(d; prob = Problem([:A, :B]), data))
end

@testset "diagnosis: a release that exists already is on the menu" begin
    # :a2 already supports the :B the user wants; what stands in the way is
    # the user's own compat on :A, and *relax your compat on A* is the fix.
    # Asking upstream for what has shipped would be asking for nothing
    data = Dict(
        :A => PkgData([:a2, :a1], Dict(:a2 => [:B], :a1 => [:B]),
                      Dict(:a2 => Dict(:B => [:b1, :b2]),
                           :a1 => Dict(:B => [:b1]))),
        :B => PkgData([:b2, :b1], DEPS_NONE, COMP_NONE),
    )
    prob = Problem([:A]; compat = Dict(:A => [:a1], :B => [:b2]))
    d = check_diagnosis(data, prob)
    c = only(d.conflicts)
    @test Set([Action(:compat, :A)]) in Set(Set(f.actions) for f in c.fixes)
    @test isempty(c.upstream)
    @test !occursin("Upstream", sprint(show, MIME("text/plain"), d))
    @test isempty(report_problems(d; prob, data))
end

@testset "diagnosis: an upstream witness settles the rest of the page" begin
    # two conflicts that share nothing: the release asked for is one conflict's
    # own, and its witness settles the other the first way that menu offers,
    # which is the convention every witness on the page follows
    data = Dict(
        :A => PkgData([:a1], Dict(:a1 => [:C]), Dict(:a1 => Dict(:C => [:c1]))),
        :C => PkgData([:c2, :c1], DEPS_NONE, COMP_NONE),
        :E => PkgData([:e1], Dict(:e1 => [:G]), Dict(:e1 => Dict(:G => [:g1]))),
        :G => PkgData([:g2, :g1], DEPS_NONE, COMP_NONE),
    )
    prob = Problem([:A, :E]; compat = Dict(:C => [:c2], :G => [:g2]))
    d = check_diagnosis(data, prob)
    @test length(d.conflicts) == 2
    for (i, c) in enumerate(d.conflicts)
        u = only(c.upstream)
        other = d.conflicts[3-i]
        # the release is this conflict's own ...
        @test u.pkg in keys(c.versions) && u.dep in keys(c.versions)
        @test u.solution[u.pkg] == u.latest
        @test u.solution[u.dep] == u.supports
        # ... and its witness is what a resolve of that registry answers with
        # the other conflict settled the first way its menu offers -- read off
        # here rather than taken from the diagnosis, and resolved from scratch
        released = Diagnostics.without_bound(data[u.pkg], u.latest, u.dep)
        mod = merge(data, Dict(u.pkg => released))
        @test u.solution == fix_resolve(mod, prob, first(other.fixes).actions)
        @test haskey(u.solution, only(other.reqs))
    end
    # ... and what prints under one is the page's own packages, not the other's
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("→ would allow: C c2\n", report)
    @test occursin("→ would allow: G g2\n", report)
    @test isempty(report_problems(d; prob, data))
end

@testset "diagnosis: the probe budget is recorded and not announced" begin
    # nine conflicts of the same shape, each with one candidate: the budget
    # stops the ninth from being tried, and what the page does not print it
    # does not talk about either
    n = Diagnostics.UPSTREAM_SOLVES + 1
    data = Dict{Symbol,PkgData{Symbol,Symbol,Vector{Symbol},Vector{Symbol},
                               Dict{Symbol,Vector{Symbol}},
                               Dict{Symbol,Dict{Symbol,Vector{Symbol}}}}}()
    reqs = Symbol[]
    compat = Dict{Symbol,Vector{Symbol}}()
    for i = 1:n
        a, b = Symbol("A", i), Symbol("B", i)
        b1, b2 = Symbol("b", i, "1"), Symbol("b", i, "2")
        data[a] = PkgData([Symbol("a", i)], Dict(Symbol("a", i) => [b]),
                          Dict(Symbol("a", i) => Dict(b => [b1])))
        data[b] = PkgData([b2, b1], DEPS_NONE, COMP_NONE)
        push!(reqs, a)
        compat[b] = [b2]
    end
    prob = Problem(reqs; compat)
    d = resolve(data, prob)
    @test length(d.conflicts) == n
    @test count(c -> !isempty(c.upstream), d.conflicts) ==
          Diagnostics.UPSTREAM_SOLVES
    @test d.upstream_cut
    report = sprint(show, MIME("text/plain"), d)
    @test !occursin("cut", report) && !occursin("budget", report)
    @test count("Upstream fix:", report) == Diagnostics.UPSTREAM_SOLVES
    @test isempty(report_problems(d; prob, data))
end

# Two qualifying pairs print as the choice they are: a bulleted list under one
# heading, each bullet the same sentence and its own witness. Built by hand,
# since a conflict that one release fixes in two different ways is rare enough
# that neither corpus holds one -- what is being checked here is the shape of
# the page, which is this file's business either way.
@testset "diagnosis: two releases print as two bullets" begin
    P, V = String, Int
    VS = Dict("X" => [1], "Q" => [1, 2, 3], "R" => [1, 2])
    lines = Line{P}[Line{P}(Clauses.clause(["X" => literal(1, [1])]), P[], true)]
    up(p, v, q, s, sup, sol) = Upstream{P,V}(p, v, q, s, V[sup...], sol)
    c = Conflict{P,V}(P["X"], lines, VS, Dict{P,Vector{Vector{Symbol}}}(),
        Fix{P,V}[Fix{P,V}(Action{P}[Action(:drop, "X")], Dict("X" => 1))],
        Tuple{Vector{Vector{Action{P}}},Vector{Action{P}}}[],
        Upstream{P,V}[up("A", 5, "Q", 3, [1, 2],
                         Dict("A" => 5, "Q" => 3, "X" => 1)),
                      up("B", 2, "R", 2, [1],
                         Dict("B" => 2, "R" => 2, "X" => 1))])
    report = sprint(show, MIME("text/plain"), Diagnosis([c], :none))
    flat = replace(report, r"\n\s+" => " ")
    @test occursin("  Upstream fixes:\n", report)
    @test !occursin("Upstream fix:", report)
    @test occursin("• a release of A supporting Q 3 would fix this; 5, its " *
                   "latest, supports only ≤2.", flat)
    @test occursin("• a release of B supporting R 2 would fix this; 2, its " *
                   "latest, supports only 1.", flat)
    # each bullet's witness is its own, and says nothing about the package the
    # sentence above it has just named a version of
    @test occursin("→ would allow: Q 3, X 1", flat)
    @test occursin("→ would allow: R 2, X 1", flat)
end
