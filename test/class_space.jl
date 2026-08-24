# The three preconditions the representation has to meet, and the primitives
# they exist for.
#
# The partition itself is tested in classes.jl, and the answers it produces are
# tested against the brute force everywhere. What is checked here is the shape
# of the representation — the properties that make a *relaxation* of a user
# constraint a well-posed question, which is why this lands before the
# diagnostics that will ask one:
#
#   1. user constraints never enter the matrices;
#   2. deactivated classes are retained as pressure-treated columns
#      (reachability continues past one, the packages behind one stay, and
#      redundancy never deletes one);
#   3. domination is an activation-independent test.
#
# Nothing here is a diagnostic. The last testsets exercise the raw primitives —
# flipping a class's activation by assumption, and probing what a reactivated
# class finds waiting for it — as operations on a built instance.

using Resolver: PkgData, PkgInfo, Problem, SAT, PicoSAT, pkg_info, nclasses,
    rank_pkg_info, prepare_pkg_info, filter_pkg_info!, deactivations,
    mark_installable!, mark_necessary!, drop_unmarked!, class_ranking,
    finalize, sat_assume, sat_pop, is_satisfiable, sat_solve

const NODEPS = Dict{Symbol,Vector{Symbol}}()
const NOCOMP = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()

# A class named by what it holds rather than by where it sits: the set of
# version *values* of its members. Class ids are per-query (the rank order is),
# so this is how two universes over one artifact are compared.
class_key(univ, p, c) =
    (p, Set(univ.info[p].versions[i] for i in univ.info[p].members[c]))

all_classes(univ) = Set{Any}(
    class_key(univ, p, c) for (p, ip) in univ.info for c = 1:nclasses(ip))

# the classes a filtering pass has struck from the universe (flag cleared)
struck(univ) = Set{Any}(
    class_key(univ, p, c)
    for (p, ip) in univ.info for c = 1:nclasses(ip) if !ip.conflicts[c, end])

# the classes this query has emptied
emptied(univ) = Set{Any}(
    class_key(univ, p, c)
    for (p, rs) in univ.reps for c in eachindex(rs) if iszero(rs[c]))

# every fact the matrices state, named so that the class *order* cannot show
# through: which classes there are, what each depends on, and which pairs
# conflict
function universe_facts(univ)
    facts = Set{Any}()
    for (p, ip) in univ.info
        for c = 1:nclasses(ip)
            key = class_key(univ, p, c)
            push!(facts, (:class, key))
            for (j, q) in enumerate(ip.depends)
                ip.conflicts[c, j] && push!(facts, (:depends, key, q))
            end
            for (q, b) in ip.interacts
                for d = 1:nclasses(univ.info[q])
                    ip.conflicts[c, b + d] || continue
                    push!(facts, (:conflict, key, class_key(univ, q, d)))
                end
            end
        end
    end
    return facts
end

# the redundancy pass on its own, on an unfiltered universe: the
# installability prune first (it is what makes the deletions safe), then
# domination
function redundancy_only(info, prob)
    univ = rank_pkg_info(info, prob)
    mark_installable!(univ.info)
    mark_necessary!(univ.info, deactivations(univ))
    return univ
end

@testset "precondition 1: no user constraint reaches a matrix" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # A constrained query and an unconstrained one make the same universe out
    # of an artifact, up to the class order the representatives put the rows
    # in: the same partition, the same rows and columns, the same conflicts.
    # Everything the constraint does is in `reps`.
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:25
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            info = pkg_info(data, keys(data); filter = false)
            base = rank_pkg_info(info, Problem(reqs))
            facts = universe_facts(base)
            @test isempty(emptied(base)) # nothing constrains, nothing is empty
            for cs = 1:4
                compat, pin = random_constraints(m, n)
                test = iseven(cs) ? (; test = (p, v) -> v == rand(1:n)) : (;)
                prob = Problem(reqs; compat, pin, test...)
                univ = rank_pkg_info(info, prob)
                # the partition is not refined: same classes, same members
                @test all_classes(univ) == all_classes(base)
                @test universe_facts(univ) == facts
                # ... and no column has appeared to carry the constraint
                @test all(size(univ.info[p].conflicts) ==
                          size(base.info[p].conflicts) for p in keys(info))
            end
        end
    end
end

@testset "precondition 2a: reachability continues past a dead class" begin
    # :P's best class is emptied by the user's compat. Stopping the prefix there
    # would keep only a class that cannot be selected, and :P would look
    # uninstallable — when in fact it resolves at its next class.
    data = Dict(
        :P => PkgData([:v2, :v1], Dict(:v2 => [:Q]), NOCOMP),
        :Q => PkgData([:w1], NODEPS, NOCOMP),
    )
    prob = Problem([:P]; compat = Dict(:P => [:v1]))
    @test resolve(data, prob) == Dict(:P => :v1)
    univ = prepare_pkg_info(pkg_info(data, prob), prob)
    # both classes survive, the first of them deactivated
    @test nclasses(univ.info[:P]) == 2
    @test univ.reps[:P] == [0, 2]

    # precondition 2b, on the same universe: :Q is depended on by nothing but
    # the deactivated class, and it is in the universe all the same — it is
    # exactly the cone a relaxation of that compat would have to find waiting
    @test haskey(univ.info, :Q)
    @test nclasses(univ.info[:Q]) == 1
end

@testset "precondition 2c: redundancy never deletes a dead class" begin
    # :P's classes are {:v3, :v1} (no constraints at all) and {:v2} (conflicts
    # with :Q@:w2), so the first dominates the second outright. Two different
    # kinds empty them: the compat bound takes the dominator, an admission
    # kind takes the dominated one. Deleting the dominated class for being
    # dominated would make relaxing the kind alone unanswerable, because the
    # class it would return is the one that was deleted.
    data = Dict(
        :P => PkgData([:v3, :v2, :v1], NODEPS, Dict(:v2 => Dict(:Q => [:w1]))),
        :Q => PkgData([:w2, :w1], NODEPS, NOCOMP),
    )
    info = pkg_info(data, [:P, :Q]; filter = false)
    @test info[:P].members == [[1, 3], [2]]
    prob = Problem([:P, :Q];
        compat = Dict(:P => [:v2]),          # empties {:v3, :v1}
        pre = (p, v) -> v == :v2)            # empties {:v2}
    univ = prepare_pkg_info(pkg_info(data, prob), prob)
    @test nclasses(univ.info[:P]) == 2
    @test univ.reps[:P] == [0, 0]
    @test emptied(univ) == Set([class_key(univ, :P, 1), class_key(univ, :P, 2)])

    sat = SAT(univ)
    try
        @test !is_satisfiable(sat, [:P])
        # relax everything, then re-deactivate only what the compat emptied:
        # the class the *kind* emptied is still there to be returned
        sat_pop(sat)
        sat_assume(sat, :P, 1; not = true)
        @test is_satisfiable(sat, [:P])
        @test PicoSAT.deref(sat.pico, sat.vars[:P] + 2) > 0
    finally
        finalize(sat)
    end
end

@testset "precondition 2c: swept" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # Two properties of the redundancy pass, over the tiny grids with random
    # constraints on top.
    #
    # First, unconditionally: it never strikes a class this query has emptied.
    # That is precondition 2's third obligation, and it holds whatever the
    # query is.
    #
    # Second, for the queries that leave the class order alone: everything the
    # pass strikes it would have struck with no constraints at all. This is the
    # general form of "no deletion rests on anything a relaxation can change",
    # and it is stated for the unreordered queries because the class *order* is
    # the one thing about a universe that a constraint legitimately does move
    # (§3: a class competes at its best admissible member, so forbidding that
    # member can put the class behind one it otherwise outranks), and domination
    # is a statement about the order — "the better class will be chosen
    # instead". A reordering query can therefore license a deletion the
    # unconstrained order does not, which is sound for the query it was
    # computed for and is what the resolve path needs; whether it is enough for
    # a later relaxation is a question about relaxation, not about this pass.
    protected = reordered = 0
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (2, 5))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:25
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            info = pkg_info(data, keys(data); filter = false)
            base = struck(redundancy_only(info, Problem(reqs)))
            for _ = 1:4
                compat, pin = random_constraints(m, n)
                prob = Problem(reqs; compat, pin)
                univ = redundancy_only(info, prob)
                gone = struck(univ)
                @test isempty(gone ∩ emptied(univ))
                protected += length(emptied(univ) ∩ base)
                if last(class_ranking(info, prob)) === nothing
                    @test gone ⊆ base
                else
                    reordered += 1
                end
            end
        end
    end
    # the sweep really did keep classes the unconstrained pass would strike ...
    @test protected > 0
    # ... and really did exercise both branches
    @test reordered > 0
end

@testset "precondition 3: domination ignores activation" begin
    # :P's best class conflicts with :Q@:w2 and its second has no constraints
    # at all, so ignoring that one conflict — and only then — would make the
    # first class's row a subset of the second's and delete the second. The
    # query empties :Q's {:w2} class, and the conflict counts all the same,
    # because what domination compares is the registry's rows and an emptied
    # class is still a class. Were it not so, relaxing the compat on :Q would
    # reactivate {:w2}, need :P@:v1, and find it deleted.
    data = Dict(
        :P => PkgData([:v2, :v1], NODEPS, Dict(:v2 => Dict(:Q => [:w1]))),
        :Q => PkgData([:w2, :w1], NODEPS, NOCOMP),
    )
    info = pkg_info(data, [:P, :Q]; filter = false)
    @test nclasses(info[:P]) == nclasses(info[:Q]) == 2

    prob = Problem([:P, :Q]; compat = Dict(:Q => [:w1]))
    univ = redundancy_only(info, prob)
    @test univ.reps[:Q] == [0, 2]
    @test isempty(struck(univ)) # :P's second class is not dominated

    # the unconstrained universe reaches the same verdict, which is the
    # property: the pass cannot tell the two queries apart
    @test struck(univ) == struck(redundancy_only(info, Problem([:P, :Q])))

    # and the end-to-end consequence: the class is still there to be reached
    # once the constraint that emptied its partner is relaxed
    full = prepare_pkg_info(pkg_info(data, prob), prob)
    sat = SAT(full)
    try
        # under the constraint, :Q is at :w1 and :P at its best
        @test is_satisfiable(sat, [:P, :Q])
        sat_pop(sat)
        # relaxed, :Q@:w2 is choosable — and :P has a class that admits it
        sat_assume(sat, :Q, 1)
        @test is_satisfiable(sat, [:P, :Q])
        @test PicoSAT.deref(sat.pico, sat.vars[:P] + 2) > 0
    finally
        finalize(sat)
    end
end

@testset "fitness: activation flips by assumption" begin
    # A class's activation literal is the negation of its class variable, so
    # flipping it is an assumption and nothing else: no clause is added, no
    # clause is rewritten, and the instance is exactly as reusable afterwards.
    data = Dict(
        :P => PkgData([:v3, :v2, :v1], NODEPS, NOCOMP),
        :Q => PkgData([:w1], Dict(:w1 => [:P]), Dict(:w1 => Dict(:P => [:v2]))),
    )
    info = pkg_info(data, [:Q])
    univ = prepare_pkg_info(info, Problem([:Q]))
    sat = SAT(univ)
    try
        clauses = PicoSAT.clause_count(sat.pico)
        @test isempty(sat.deact) # nothing is deactivated to begin with
        @test is_satisfiable(sat, [:Q])
        # :Q needs :P at the class holding :v2; deactivating that class by
        # assumption makes the requirement unsatisfiable ...
        c = univ.info[:P].classes[findfirst(==(:v2), univ.info[:P].versions)]
        sat_assume(sat, :P, c; not = true)
        @test !is_satisfiable(sat, [:Q])
        # ... and the very next solve sees it active again, on an instance with
        # exactly the clauses it started with
        @test is_satisfiable(sat, [:Q])
        @test PicoSAT.clause_count(sat.pico) == clauses
    finally
        finalize(sat)
    end
end

@testset "fitness: the D1′ shape is unrepresentable" begin
    # Proposition D1′'s counterexample shape: a compat bound and a pin that
    # between them forbid both versions of a package, where relaxing one of the
    # two kinds ought to bring one version back. What made a partial
    # relaxation unsound in version space was that the constraints entered the
    # rows, could make the two versions row-identical, and redundancy
    # elimination could then delete the one the relaxation needed.
    #
    # Here the constraints are not rows, and versions with identical registry
    # rows are one class sharing one column. There is no pair of separately
    # deletable objects for a partial relaxation to fall between.
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]), NOCOMP),
        :B => PkgData([:w1], NODEPS, NOCOMP),
    )
    prob = Problem([:A]; compat = Dict(:A => [:v1]), pin = Dict(:A => :v2))
    info = pkg_info(data, prob)
    # nothing in the registry tells :v2 and :v1 apart, so they are one class
    @test info[:A].classes == [1, 1]
    @test info[:A].members == [[1, 2]]
    @test nclasses(info[:A]) == 1
    # the two kinds empty that one class between them
    univ = rank_pkg_info(info, prob)
    @test univ.reps[:A] == [0]
    @test resolve(data, prob; diagnose = false) === nothing

    # relaxing either kind reactivates the same single class — there is no
    # "half" of this constraint set to be unsound about, and nothing that could
    # have been deleted for being dominated by the other half
    for (relaxed, want) in ((Problem([:A]; pin = Dict(:A => :v2)), :v2),
                            (Problem([:A]; compat = Dict(:A => [:v1])), :v1))
        u = rank_pkg_info(info, relaxed)
        @test nclasses(u.info[:A]) == 1
        @test u.info[:A].members == [[1, 2]]
        @test u.reps[:A] != [0]
        @test resolve(data, relaxed) == Dict(:A => want, :B => :w1)
    end

    # the same in the SAT instance: one variable, one deactivating clause, and
    # popping the frame is the whole relaxation
    sat = SAT(univ)
    try
        @test sat.deact == [-(sat.vars[:A] + 1)]
        @test !is_satisfiable(sat, [:A])
        sat_pop(sat)
        @test is_satisfiable(sat, [:A])
    finally
        finalize(sat)
    end
end

@testset "fitness: a reactivated class finds its cone" begin
    # Precondition 2b end to end. :P@:v2 is the only thing in the registry that
    # depends on :R, and the query's compat empties the class holding it. The
    # universe still has to carry :R — otherwise "what if that compat were
    # relaxed?" reactivates a class into a universe with nothing to satisfy it.
    data = Dict(
        :P => PkgData([:v2, :v1], Dict(:v2 => [:R]), NOCOMP),
        :R => PkgData([:s2, :s1], NODEPS, NOCOMP),
    )
    prob = Problem([:P]; compat = Dict(:P => [:v1]))
    univ = prepare_pkg_info(pkg_info(data, prob), prob)
    @test univ.reps[:P] == [0, 2]
    @test haskey(univ.info, :R)

    sat = SAT(univ)
    try
        # :R has a variable of its own, and the dependency clause is there
        @test haskey(sat.vars, :R)
        @test is_satisfiable(sat, [:P])
        # ... and :P's first class, the only thing that needs it, is out
        sat_assume(sat, :P, 1)
        @test !sat_solve(sat)
        # reactivate, and the cone is right there: choosing :P's first class
        # forces :R in
        sat_pop(sat)
        sat_assume(sat, :P, 1)
        @test sat_solve(sat)
        @test PicoSAT.deref(sat.pico, sat.vars[:R]) > 0
    finally
        finalize(sat)
    end
end

# the versions a universe can still name but no longer offers
shadowed(univ, p) = Set(v for sh in univ.info[p].shadows for v in sh)

# the redundancy pass on its own, and the classes it struck: everything the
# installability prune took is subtracted, so what is left is exactly the
# deletions that owe an explanation
function redundancy_delta(info, prob)
    univ = rank_pkg_info(info, prob)
    mark_installable!(univ.info)
    before = struck(univ)
    mark_necessary!(univ.info, deactivations(univ), univ.ranks)
    return univ, setdiff(struck(univ), before)
end

@testset "shadows: a deleted version keeps a place to be named from" begin
    # :P's four versions are four classes — {:Q}, {:Q,:R}, {:S}, {:S,:R} — and
    # two of them dominate: :v4's row is a subset of :v3's, :v2's of :v1's. The
    # two survivors are not comparable to each other, so each has to keep its
    # own deletion rather than the universe keeping one pile of them
    data = Dict(
        :P => PkgData([:v4, :v3, :v2, :v1],
            Dict(:v4 => [:Q], :v3 => [:Q, :R],
                 :v2 => [:S], :v1 => [:S, :R]), NOCOMP),
        :Q => PkgData([:w1], NODEPS, NOCOMP),
        :R => PkgData([:w1], NODEPS, NOCOMP),
        :S => PkgData([:w1], NODEPS, NOCOMP),
    )
    info = pkg_info(data, [:P]; filter = false)
    @test nclasses(info[:P]) == 4
    univ, _ = redundancy_delta(info, Problem([:P]))
    drop_unmarked!(univ)

    # what survives, and what each survivor speaks for
    @test univ.info[:P].versions == [:v4, :v2]
    @test univ.info[:P].members == [[1], [2]]
    @test univ.info[:P].shadows == [[:v3], [:v1]]
    # every version is accounted for: offered or shadowed, never lost
    @test Set(univ.info[:P].versions) ∪ shadowed(univ, :P) ==
          Set(data[:P].versions)

    # the whole filter then walks past :v2's class — nothing requires :S — and
    # the shadow goes with it, which is what shadowing a *class* means
    full = prepare_pkg_info(info, Problem([:P]))
    @test full.info[:P].versions == [:v4]
    @test full.info[:P].shadows == [[:v3]]
    # the artifact is left alone: a universe shares its shadow lists with the
    # artifact it was copied from, so growing one has to replace it, and a
    # second resolve over the same artifact has to see the same universe
    @test prepare_pkg_info(info, Problem([:P])).info == full.info
    # ... and the answer is the one it was before any of this
    @test resolve(data, [:P]) == Dict(:P => :v4, :Q => :w1)

    # a chain lands on the class that survives it, not on the nearest one:
    # :v3 dominates :v2, but :v4 dominates both
    data = Dict(
        :P => PkgData([:v4, :v3, :v2, :v1],
            Dict(:v3 => [:Q], :v2 => [:Q, :R], :v1 => [:Q, :R, :S]), NOCOMP),
        :Q => PkgData([:w1], NODEPS, NOCOMP),
        :R => PkgData([:w1], NODEPS, NOCOMP),
        :S => PkgData([:w1], NODEPS, NOCOMP),
    )
    univ = prepare_pkg_info(pkg_info(data, [:P]; filter = false),
                            Problem([:P]))
    @test univ.info[:P].versions == [:v4]
    @test univ.info[:P].shadows == [[:v3, :v2, :v1]]
    @test resolve(data, [:P]) == Dict(:P => :v4)
end

@testset "shadows: a shadow is not a member" begin
    # :v1 needs :Q and :v2 does not, so :v2's class dominates :v1's and, with
    # nothing constraining :P, :v1 is deleted and shadowed by it.
    data = Dict(
        :P => PkgData([:v2, :v1], Dict(:v1 => [:Q]), NOCOMP),
        :Q => PkgData([:w1], NODEPS, NOCOMP),
    )
    info = pkg_info(data, [:P]; filter = false)
    free = prepare_pkg_info(info, Problem([:P]))
    @test free.info[:P].versions == [:v2]
    @test free.info[:P].shadows == [[:v1]]
    @test !haskey(free.info, :Q) # nothing surviving pulls it in

    # Now pin :P to the shadowed version. Were the shadow a *member* of the
    # class that shadows it, the class would simply stand at :v1 — and the
    # class's row says nothing about :Q, so the answer would leave :Q out and
    # be wrong. A shadow is not a member: the pin empties the dominating class,
    # an emptied class dominates nothing, and :v1's own class — with its own
    # row, and its own dependency — is there to be selected.
    prob = Problem([:P]; pin = Dict(:P => :v1))
    univ = prepare_pkg_info(info, prob)
    @test univ.info[:P].versions == [:v2, :v1]
    @test univ.info[:P].members == [[1], [2]]
    @test all(isempty, univ.info[:P].shadows)
    @test univ.reps[:P] == [0, 2] # the dominator is the emptied one
    @test resolve(data, prob) == Dict(:P => :v1, :Q => :w1)
    test_bake_equivalence(data, prob)
end

@testset "shadows: swept" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # Two properties of the shadow lists over the tiny grids with random
    # constraints on top.
    #
    # First, they account for the redundancy pass exactly: every version that
    # pass deleted is named by some surviving class, and no other version is.
    # The grids give the installability prune nothing to do — every package has
    # versions and every dependency names one of them — which is what makes the
    # accounting exact rather than an inclusion, so that is asserted too.
    #
    # Second, the claim the lists exist to support: the class a version is
    # shadowed by has a *subset* of its constraints. That is what makes "this
    # class is ruled out" carry over to everything it shadows, and it is
    # checked here against the raw data rather than against the matrix the
    # pass read, so a bug in the row test cannot hide in both.
    found = 0
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (2, 5))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:25
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            info = pkg_info(data, keys(data); filter = false)
            for _ = 1:4
                compat, pin = random_constraints(m, n)
                prob = Problem(reqs; compat, pin)
                univ, delta = redundancy_delta(info, prob)
                pruned = Dict(p => Set(ip.versions[i]
                        for cl = 1:nclasses(ip) if !ip.conflicts[cl, end]
                        for i in ip.members[cl])
                    for (p, ip) in univ.info)
                gone = Dict{Any,Set}()
                for (p, vs) in delta
                    union!(get!(() -> Set(), gone, p), vs)
                    setdiff!(pruned[p], vs)
                end
                drop_unmarked!(univ)
                for (p, ip) in univ.info
                    sh = shadowed(univ, p)
                    @test isempty(pruned[p])
                    @test sh == get(() -> Set(), gone, p)
                    found += length(sh)
                    # the subset claim, read off the raw data
                    for cl = 1:nclasses(ip), v in ip.shadows[cl]
                        r = ip.versions[first(ip.members[cl])]
                        @test Set(lookup(data[p].depends, r, ())) ⊆
                              Set(lookup(data[p].depends, v, ()))
                        for (q, iq) in univ.info, w in iq.versions
                            forbids(u) = ref_forbids(data, p, u, q, w) ||
                                         ref_forbids(data, q, w, p, u)
                            forbids(r) && @test forbids(v)
                        end
                    end
                end
            end
        end
    end
    @test found > 0 # the sweep really did find some
end
