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
using Resolver.Diagnostics: Diagnostics, Conflict, Fix, Action, Fact,
    Requirement, Availability, Dependency, Incompatibility, unavailable,
    action_phrase
using Resolver.UnsatCores: sat_mcses

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

# A dependency and an incompatibility claim none: they say what the registry
# says, which is not something this query could relax, so there is nothing for
# the solver to be asked about them and nothing taking one away would give
# back. What they say is checked against the universe instead, by
# `relations_are_true` below
fact_literals(::SAT{P}, ::Dependency{P}) where {P} = Int[]
fact_literals(::SAT{P}, ::Incompatibility{P}) where {P} = Int[]

conflict_literals(sat::SAT, c::Conflict) =
    reduce(vcat, (fact_literals(sat, f) for f in c.chain); init = Int[])

# is the instance satisfiable assuming exactly `lits`?
function sat_assuming(sat::SAT, lits)
    for l in lits
        sat_assume_var(sat, l)
    end
    return sat_solve(sat)
end

# Does the chain consist of exactly the reasons its requirements fail? Every
# subset of its facts is put to the solver, the minimal unsatisfiable ones are
# read off, and what has to hold is that between them they use every fact.
#
# This is the general form of "every fact is load-bearing", which is the case
# where the only minimal subset is the whole chain — a conflict whose
# requirements fail *together*. A conflict that gathers requirements failing
# separately for a shared reason has one minimal subset per requirement, and
# the test is that the chain is their union: no fact left over, and none of the
# requirements dropped for being redundant.
#
# Brute force, which is what makes it an oracle; chains are a handful of facts,
# and one too big to enumerate falls back to the unsatisfiability check alone.
function chain_is_covered(sat::SAT, c::Conflict)
    # only the facts that claim literals are enumerated over: the rest say what
    # the registry says, which no subset of assumptions can turn off, so they
    # are in every subset and in none of them alike
    lits = [fact_literals(sat, f) for f in c.chain]
    claim = Int[i for i in eachindex(lits) if !isempty(lits[i])]
    n = length(claim)
    n ≤ 8 || return !sat_assuming(sat, conflict_literals(sat, c))
    subset(bits) = reduce(vcat,
        (lits[claim[i]] for i = 1:n if !iszero(bits & (1 << (i-1))));
        init = Int[])
    unsat = Bool[!sat_assuming(sat, subset(bits)) for bits = 0:(1 << n) - 1]
    covered = falses(n)
    for bits = 0:(1 << n) - 1
        unsat[bits+1] || continue
        # minimal iff dropping any one of its facts makes it satisfiable
        any(i -> !iszero(bits & (1 << (i-1))) && unsat[(bits ⊻ (1 << (i-1)))+1],
            1:n) && continue
        for i = 1:n
            iszero(bits & (1 << (i-1))) || (covered[i] = true)
        end
    end
    return all(covered)
end

# What one class of `p` says about `q`: whether it depends on it, and which of
# `q`'s versions it admits. Read straight off the universe's matrices, which is
# where a relation fact's claim has to be checked — the solver is never asked
# about a dependency edge or a registry bound, so it has nothing to say here
function class_relation(sat::SAT{P}, p::P, c::Int, q::P) where {P}
    info_p, info_q = sat.info[p], sat.info[q]
    k = findfirst(==(q), info_p.depends)
    dep = k !== nothing && info_p.conflicts[c, k]
    allowed = fill(true, length(info_q.versions))
    if haskey(info_p.interacts, q)
        b = info_p.interacts[q]
        for j = 1:nclasses(info_q), i in info_q.members[j]
            allowed[i] = !info_p.conflicts[c, b + j]
        end
    end
    return dep, allowed
end

# the versions of `p` a fact speaks for, as indices into what `p` offers
fact_indices(sat::SAT{P}, p::P, versions) where {P} =
    Int[findfirst(==(v), sat.info[p].versions) for v in versions]

# Does every relation in the chain say something true of the universe, in the
# only direction it may be said in?
#
#   * a dependency is stated only where the registry has one. `PkgInfo` records
#     compatibility symmetrically, so a bound cannot be attributed to a package
#     that merely has one — that is what the other kind is for;
#   * an incompatibility is stated where the versions it speaks for do not
#     depend on what it speaks about, so there is no side to put the bound on
#     and it puts it on neither. Reading it the other way round says the same
#     thing, since that is what symmetric means;
#   * the versions it speaks for are versions the query admits, they are a run
#     of what the package offers, and every one of them says exactly what the
#     fact says. A fact that averaged over versions that disagree would be a
#     bound none of them declares.
function relations_are_true(sat::SAT{P}, prob::Problem{P}, c::Conflict) where {P}
    for f in c.chain
        f isa Dependency || f isa Incompatibility || continue
        p = f.pkg
        q = f isa Dependency ? f.dep : f.other
        info_p, info_q = sat.info[p], sat.info[q]
        @test f.offering == info_q.versions
        idx = fact_indices(sat, p, f.versions)
        @test idx == first(idx):last(idx)  # a run of what `p` offers
        for i in idx
            v = info_p.versions[i]
            @test isempty(exclusion_kinds(prob, p, v))
            dep, allowed = class_relation(sat, p, info_p.classes[i], q)
            @test dep == (f isa Dependency)
            f.allowed === nothing || @test allowed == f.allowed
        end
        # the other way round: reading an incompatibility from `q`'s side picks
        # out the same pairs
        f isa Incompatibility || continue
        for (j, u) in enumerate(info_q.versions)
            _, back = class_relation(sat, q, info_q.classes[j], p)
            @test f.allowed[j] == all(back[i] for i in idx)
        end
    end
end

# A relation with no bound says only that the package is needed. That is
# enough where something else in the story says what would have done instead —
# the next step, whose subject is the versions this one would have named, or
# the relation the story ends on — and enough where the query has left the
# package nothing at all, since then every bound would be starved alike. It is
# not enough anywhere else: dropping a bound there would drop the answer
function bounds_are_stated(sat::SAT{P}, c::Conflict{P}) where {P}
    rels = Fact[f for f in c.chain
                if f isa Dependency || f isa Incompatibility]
    # a package the story goes on from, or that some bound speaks about
    spoken = Set{P}(f.pkg for f in rels)
    for f in rels
        f.allowed === nothing ||
            push!(spoken, f isa Dependency ? f.dep : f.other)
    end
    return all(rels) do f
        f isa Dependency && f.allowed === nothing || return true
        return f.dep in spoken || all(iszero, sat.reps[f.dep])
    end
end

# Do the relations in a chain hang together — every package one speaks of
# either something the conflict requires, or something an earlier relation
# already brought in? A relation about packages nothing else in the story
# mentions would be true and beside the point
function relations_connect(c::Conflict{P}) where {P}
    reached = Set{P}(c.reqs)
    for f in c.chain
        f isa Dependency || f isa Incompatibility || continue
        f.pkg in reached || return false
        push!(reached, f isa Dependency ? f.dep : f.other)
    end
    return true
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

        # each names requirements of the query, and its chain names them in
        # order and then the packages
        for c in d.conflicts
            @test !isempty(c.reqs)
            @test c.reqs ⊆ prob.reqs
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
                # ... and nothing in it is a passenger: the minimal conflicts
                # among its own facts cover every one of them
                @test chain_is_covered(sat, c)
            end
        end

        # the registry's own facts say what the registry says, in a direction
        # it licenses, about packages the rest of the story reaches
        for c in d.conflicts
            relations_are_true(sat, prob, c)
            @test relations_connect(c)
            @test bounds_are_stated(sat, c)
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

@testset "diagnosis: the story names the constraint" begin
    prob = Problem([:A]; compat = Dict(:B => [:w1]), pin = Dict(:B => :w2))
    d = check_diagnosis(needs_dep, prob)
    @test length(d.conflicts) == 1
    c = only(d.conflicts)
    @test c.reqs == [:A]
    # the middle fact is the one the user did not write: :A needs :B, and no
    # bound of :A's is what the story turns on, since the query has left :B
    # with nothing whatever :A would have taken
    @test c.chain == Fact[Requirement(:A),
        Dependency{Symbol,Symbol}(:A, [:v1], :B, [:w3, :w2, :w1], nothing, true, true),
        Availability{Symbol,Symbol}(:B, [:w3, :w2, :w1],
            [[:compat, :pin], [:compat], [:pin]])]
    @test unavailable(c.chain[3])
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
    @test d.conflicts[1].chain[2] ==
        Dependency{Symbol,Symbol}(:R, [:r1], :P, [:p2, :p1], [true, false], true, true)
    f = d.conflicts[1].chain[3]::Availability
    @test f.pkg == :P
    @test f.members == [:p2, :p1]
    @test f.excluded == [[:compat], Symbol[]]
    @test !unavailable(f)
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict:

        Conflict 1: R cannot be satisfied.
          • you require R
          • R r1 requires P at p2
          • P: p2 excluded by your compat
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
    @test c.chain == Fact[Requirement(:A),
        Dependency{Symbol,Symbol}(:A, [:a3], :C, [:c3, :c2, :c1],
            [true, false, false], true, false),
        Dependency{Symbol,Symbol}(:A, [:a2, :a1], :C, [:c3, :c2, :c1],
            [false, true, false], false, true),
        Availability{Symbol,Symbol}(:C, [:c3, :c2, :c1],
            [[:compat], [:compat], Symbol[]])]
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict:

        Conflict 1: A cannot be satisfied.
          • you require A
          • A a3 requires C at c3
          • A ≤ a2 requires C at c2
          • C: ≥ c2 excluded by your compat
          Fix it by any one of:
            1. relax your compat on C
               → allows: A a3, C c3
            2. drop requirement A
        """
end

@testset "diagnosis: a story that spans more than one hop" begin
    # the query names :A and :C, and it is :B in between that explains them
    d = check_diagnosis(two_hops, Problem([:A]; compat = Dict(:C => [:c1])))
    c = only(d.conflicts)
    # the step to :B says only that :B is needed. What :A will take of it is
    # not what the query left nothing of, and the versions it would name are
    # the subject of the next step anyway
    @test c.chain == Fact[Requirement(:A),
        Dependency{Symbol,Symbol}(:A, [:a1], :B, [:b1], nothing, true, true),
        Dependency{Symbol,Symbol}(:B, [:b1], :C, [:c2, :c1], [true, false], true, true),
        Availability{Symbol,Symbol}(:C, [:c2, :c1], [[:compat], Symbol[]])]
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("A a1 requires B\n", report)
    @test occursin("B b1 requires C at c2\n", report)
    # ... and :B is not something the query said anything about: the middle of
    # the story is the part only the registry knows
    @test :B ∉ [f.pkg for f in c.chain
                if f isa Requirement || f isa Availability]
end

@testset "diagnosis: a bound with no dependency behind it" begin
    # :P's bound on :W is a weak dependency: a compat entry with no edge. The
    # universe records compatibility symmetrically, so which of the two
    # declared it is not recoverable — and saying ":P requires :W" would be
    # attributing a bound to a package that does not even depend on the other
    prob = Problem([:P, :S]; compat = Dict(:W => [:w2]))
    d = check_diagnosis(weak_bound, prob)
    c = only(d.conflicts)
    @test c.chain == Fact[Requirement(:P), Requirement(:S),
        Dependency{Symbol,Symbol}(:S, [:s1], :W, [:w2, :w1], nothing, true, true),
        Incompatibility{Symbol,Symbol}(:P, [:p1], :W, [:w2, :w1],
            [false, true], true, true),
        Availability{Symbol,Symbol}(:W, [:w2, :w1], [Symbol[], [:compat]])]
    # the only dependency stated is the one the registry has
    @test [(f.pkg, f.dep) for f in c.chain if f isa Dependency] == [(:S, :W)]
    @test !any(f isa Dependency && (f.pkg == :P || f.dep == :P) for f in c.chain)
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("S s1 requires W\n", report)
    @test occursin("P p1 works with W only at w1\n", report)
    @test !occursin("P p1 requires", report)
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
    f = d.conflicts[1].chain[3]::Availability
    @test f.excluded == [[:prerelease], [:prerelease], [:yanked]]
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
    @test d.conflicts[1].chain == Fact[Requirement(:A),
        Availability{Symbol,Symbol}(:A, Symbol[], Vector{Symbol}[])]
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
    @test occursin("The only fix: relax your compat on C", report)
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
    # needs :C, and the bound on :C they share
    @test c.reqs == [:A, :B]
    @test c.chain == Fact[Requirement(:A), Requirement(:B),
        Availability{Symbol,Symbol}(:A, [:a2, :a1], [Symbol[], [:compat]]),
        Availability{Symbol,Symbol}(:B, [:b2, :b1], [Symbol[], [:compat]]),
        Dependency{Symbol,Symbol}(:A, [:a2], :C, [:c2, :c1], [true, false], true, true),
        Dependency{Symbol,Symbol}(:B, [:b2], :C, [:c2, :c1], [true, false], true, true),
        Availability{Symbol,Symbol}(:C, [:c2, :c1], [[:compat], Symbol[]])]
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
    @test occursin("you require A", report)
    @test occursin("you require B", report)
    @test occursin("A and B cannot both be satisfied", report)
    # both of them need :C at the version the bound took, and naming one of
    # them would leave the other's failure unexplained
    @test occursin("A a2 requires C at c2", report)
    @test occursin("B b2 requires C at c2", report)
    @test occursin("The only fix: relax your compat on C", report)
    # the witness names both of them too
    @test occursin("→ allows: A a2, B b2, C c2", report)
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
          • A v1 requires C at v1
          • B v1 requires C at v2
          Fix it by any one of:
            1. drop requirement A
               → allows: B v1, C v2
            2. drop requirement B
               → allows: A v1, C v1

        Conflict 2: E and F cannot both be satisfied.
          • you require E
          • you require F
          • E v1 requires G at v1
          • F v1 requires G at v2
          Fix it by any one of:
            1. drop requirement E
               → allows: F v1, G v2
            2. drop requirement F
               → allows: E v1, G v1
        """
    # the one-line summary counts the ways of repairing the whole query
    @test sprint(show, d) == "Diagnosis: 2 conflicts, 4 fixes"

    # every kind that excludes a version is named, and a long line is filled
    d = resolve(needs_dep,
        Problem([:A]; compat = Dict(:B => [:w1]), pin = Dict(:B => :w2)))
    @test sprint(show, MIME("text/plain"), d) == """
        Unsatisfiable — 1 conflict:

        Conflict 1: A cannot be satisfied.
          • you require A
          • A v1 requires B
          • no version of B is available: w3 excluded by your compat and your pin, w2
            by your compat, w1 by your pin
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
    @test any(f isa Availability && f.pkg == "PrettyTables" && unavailable(f)
              for f in c.chain)
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
    f = only(x for x in c.chain if x isa Dependency)
    @test f.pkg == "DataFrames" && f.dep == "PrettyTables"
    @test f.allowed === nothing
    @test all(v -> v ∈ VersionSpec("1"), f.versions)
    report = sprint(show, MIME("text/plain"), d)
    @test occursin("no version of PrettyTables is available", report)
    @test occursin("relax your compat on PrettyTables", report)
    @test occursin(
        "DataFrames $(f.versions[end])–$(f.versions[1]) requires PrettyTables\n",
        report)

    # ... and a query that leaves PrettyTables something DataFrames will not
    # take is where the bound itself is the story
    prob = Problem(["DataFrames"];
        compat = Dict("DataFrames" => VersionSpec("1.4 - 1.7"),
                      "PrettyTables" => VersionSpec("1")))
    d = check_diagnosis(rp, prob)
    @test d isa Diagnosis{String,VersionNumber}
    f = only(x for x in only(d.conflicts).chain if x isa Dependency)
    @test f.pkg == "DataFrames" && f.dep == "PrettyTables"
    @test all(v -> v ∈ VersionSpec("1.4 - 1.7"), f.versions)
    # what it will take is PrettyTables 2, which is what the bound rules out
    @test all(v ∈ VersionSpec("2") for (v, ok) in zip(f.offering, f.allowed) if ok)
    @test any(f.allowed)
    @test occursin("requires PrettyTables at ",
        sprint(show, MIME("text/plain"), d))

    # a bound that differs across the depending package's versions, on real
    # data: Plots has wanted a different RecipesBase over the years, and a
    # query holding it at 1 and RecipesBase at 0.4 is stopped by all of them
    prob = Problem(["Plots"];
        compat = Dict("Plots" => VersionSpec("1"),
                      "RecipesBase" => VersionSpec("0.4")))
    d = check_diagnosis(rp, prob)
    @test d isa Diagnosis{String,VersionNumber}
    deps = Dependency{String,VersionNumber}[f for f in only(d.conflicts).chain
                                            if f isa Dependency]
    @test all(f -> f.pkg == "Plots" && f.dep == "RecipesBase", deps)
    # several of them, because the versions of Plots disagree — and what they
    # disagree about is exactly what makes them several
    @test length(deps) > 1
    @test allunique(f.allowed for f in deps)
    # each speaks for its own versions, and between them they run newest to
    # oldest over the versions of Plots the query admits
    spoken = vcat((f.versions for f in deps)...)
    @test allunique(spoken)
    @test issorted(spoken; rev = true)
    @test all(v -> v ∈ VersionSpec("1"), spoken)
end
