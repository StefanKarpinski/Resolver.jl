# The property the filter buys with the two rank keys: a universe filtered for
# query `Q` still answers every *relaxation* of `Q` — the same query with some
# requirements or constraint sources dropped, which is exactly the shape of a
# repair a user can make to an unsatisfiable problem.
#
# The theory page states this as Theorem C. What is checked here is its
# conclusion, end to end and without reference to how it is proved:
#
#     re-rank the Q-filtered universe under R, solve it, and you get what
#     resolving R from the unfiltered artifact gives
#
# for random universes crossed with random queries and random relaxations of
# them. The left-hand side is `resolve_relaxed`, run the way a diagnosis runs
# it: one SAT instance per query, answering the whole sample of relaxations.
#
# Two sharper statements are checked alongside, because a failure in either
# localizes the cause: that the Q-filtered universe contains everything
# a fresh filter for R would keep (Theorem B), and that the reachability walk's
# kept set is downward closed in `key_∅` (the Remark), which is the invariant
# that lets it carry an integer frontier instead of a set.

using Resolver: PkgData, Problem, SAT, NoKinds, ClassRanking,
    pkg_info, nclasses, prepare_pkg_info, rank_pkg_info, resolve,
    resolve_relaxed, class_ranking, version_permutations, deactivations,
    is_excluded, with_classes_relaxed, with_temp_clauses, sat_push,
    mark_installable!, mark_reachable!, finalize

using .tiny_data

# --- relaxations, at the Problem level: exactly what a diagnosis would ask ---

# every constraint source Q carries, as something removable
function sources(Q::Problem{P}) where {P}
    s = Any[]
    for p in keys(Q.compat); push!(s, (:compat, p)); end
    for p in keys(Q.pins);   push!(s, (:pin, p));    end
    for (k, _) in Q.excludes; push!(s, (:kind, k));  end
    return s
end

# What being a relaxation *is*, checked on every one built here. `R ≤ Q` means
# `reqs_R ⊆ reqs_Q` and `adm_Q ⊆ adm_R` — every version `Q` admitted, `R` admits
# — which is the hypothesis (Lemma 1) that every theorem on the
# relaxation-stability page rests on. `resolve_relaxed` takes it as a
# precondition and answers wrongly rather than erroring when it fails, so what
# has to be right is the construction, and this is where the construction is.
function is_relaxation(info, Q::Problem{P}, R::Problem{P}) where {P}
    Set(R.reqs) ⊆ Set(Q.reqs) || return false
    for (p, info_p) in info, v in info_p.versions
        is_excluded(R, p, v) && !is_excluded(Q, p, v) && return false
    end
    return true
end

function relax(info, Q::Problem{P}, drop_reqs, drop_srcs) where {P}
    reqs = P[p for p in Q.reqs if p ∉ drop_reqs]
    nocomp = Set(p for (k, p) in drop_srcs if k === :compat)
    nopin  = Set(p for (k, p) in drop_srcs if k === :pin)
    nokind = Set(k for (t, k) in drop_srcs if t === :kind)
    compat = filter(kv -> first(kv) ∉ nocomp, Dict(Q.compat))
    pins   = filter(kv -> first(kv) ∉ nopin,  Dict(Q.pins))
    excludes = Pair{Symbol,Any}[kv for kv in Q.excludes if first(kv) ∉ nokind]
    R = Problem(reqs; compat, pins, excludes)
    @test is_relaxation(info, Q, R)
    return R
end

# a sample: the bottom of the lattice (everything dropped) first, since it is
# the hardest case and the one caching depends on, then random points
function relaxations(info, Q::Problem{P}, k::Int) where {P}
    srcs = sources(Q)
    out = [relax(info, Q, P[], srcs)]
    for _ = 1:k
        dr = P[p for p in unique(Q.reqs) if rand() < 0.35]
        ds = Any[s for s in srcs if rand() < 0.5]
        isempty(dr) && isempty(ds) && continue
        push!(out, relax(info, Q, dr, ds))
    end
    return out
end

# Constraints that demote every class that can be demoted: forbid the best
# member of each multi-member class, so the class competes at a worse member and
# can fall behind one it used to outrank — and nothing is emptied, so this
# reorders without deactivating. Dropping the bound puts the classes back, which
# is a relaxation whose whole content is the class order moving.
function demoting_constraints(info::AbstractDict{P}) where {P}
    compat = Dict{P,Vector{Int}}()
    for (p, info_p) in info
        banned = Set(info_p.versions[first(mem)]
                     for mem in info_p.members if length(mem) ≥ 2)
        isempty(banned) && continue
        compat[p] = [v for v in info_p.versions if v ∉ banned]
    end
    return compat
end

# One round's data and constraints, in one of two regimes. `:random` is the
# ordinary case. `:demoting` builds the constraints above, and needs classes of
# several members for them to have anything to demote — which sparse dependency
# and compat patterns produce (the same trick `test/classes.jl` uses, for the
# same reason), so it draws sparsely and redraws until the data has one.
#
# Only a class with two members can move: a singleton competes at its only
# member whatever the query says, so its key is the same under every relaxation.
# That is why reordering is rare under random constraints, and why testing it
# takes a regime of its own rather than a bigger sample.
function draw_query!(mode::Symbol, m, n, make_deps, make_comp, d, c, data)
    for _ = 1:25
        sparse = mode === :demoting
        fill_data!(m, n,
            make_deps(sparse ? randbits(d) & randbits(d) : randbits(d)),
            make_comp(sparse ? randbits(c) & randbits(c) & randbits(c) :
                               randbits(c)), data)
        info = pkg_info(data, keys(data); filter = false)
        sparse || return info, random_constraints(m, n)...
        compat = demoting_constraints(info)
        isempty(compat) || return info, compat, Dict{Int,Int}()
    end
    return nothing # no data with a class of several members: nothing to demote
end

# the (package, version-tuple) identity of every class left in a universe, so
# that two independently built universes can be compared by content
classkeys(univ) = Set(
    (p, Tuple(ip.versions[j] for j in mem))
    for (p, ip) in univ.info for mem in ip.members)

# does R move any of the universe's classes, or only re-represent them? the
# first is what makes the descent optimize in an order that is not the layout
reranks(univ, R) =
    class_ranking(univ.info, R, version_permutations(univ.info, nothing))[2] !== nothing

# classes the query emptied that R gives a representative back to
revivals(univ, R) = count(
    iszero(univ.reps[p][c]) && !iszero(rs[c])
    for (p, rs) in class_ranking(univ.info, R)[1] for c in eachindex(rs))

@testset "a Q-filtered universe answers every relaxation of Q" begin
    reused = revived = reranked = demoting = 0
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (2, 5))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        # both regimes: constraints drawn at random, and constraints built to
        # demote every class that can be demoted
        for mode in (:random, :demoting), _ = 1:15
            drawn = draw_query!(mode, m, n, make_deps, make_comp, d, c, data)
            drawn === nothing && continue
            info, compat, pins = drawn
            reqs = collect(make_reqs(rand(1:2^m-1)))
            bad = rand(1:n)
            # a second source, independently droppable, so that dropping the
            # compat bound is one point of the lattice rather than all of it
            excludes = mode === :demoting || rand(Bool) ?
                Pair{Symbol,Any}[:test => (p, v) -> v == bad] : NoKinds
            Q = Problem(reqs; compat, pins, excludes)
            univ = prepare_pkg_info(info, Q)
            have = classkeys(univ)
            # one instance for the whole lattice sample, which is the point:
            # a diagnosis tests its candidate fixes on the instance the failed
            # resolve already built
            sat = SAT(univ)
            try
                for R in relaxations(info, Q, 4)
                    # Theorem C: reuse answers R exactly as a fresh resolve does
                    @test resolve_relaxed(sat, R) == resolve(info, R)
                    # Theorem B, the sharper statement: nothing a fresh filter
                    # for R keeps is missing from the universe filtered for Q
                    @test isempty(setdiff(classkeys(prepare_pkg_info(info, R)), have))
                    reused += 1
                    revived += revivals(univ, R) > 0
                    reranked += reranks(univ, R)
                    demoting += mode === :demoting
                end
            finally
                finalize(sat)
            end
        end
    end
    @test reused > 800
    # the sweep really did exercise both things a relaxation moves: the classes
    # the descent can select at all, and the order it optimizes them in
    @test revived > 0
    # the second is what the demoting regime is for. It reorders 11-15% of the
    # relaxations it sweeps, against 0.5-1% under random constraints alone —
    # dozens of times a run rather than a handful, so this is a floor with room
    # under it rather than a coin flip. It is not higher because reordering also
    # needs the demoted class to be non-contiguous in version order, so that
    # something it outranked sits between two of its members
    @test reranked > demoting ÷ 15
end

@testset "a relaxation the descent optimizes in a new order" begin
    # :A's :v4 and :v1 are one class (both depend on :B alone), :v3 and :v2 are
    # another (both on :C). Q forbids :v4, which leaves the first class standing
    # at :v1 — behind the second, which stands at :v3 — so the universe is laid
    # out with the :C class first. R lifts the bound and the order flips back:
    # "some better class" then means a class *later* in the layout
    data = Dict(
        :A => PkgData([:v4, :v3, :v2, :v1],
            Dict(:v4 => [:B], :v3 => [:C], :v2 => [:C], :v1 => [:B]), NO_COMP),
        :B => PkgData([:w1], NO_DEPS, NO_COMP),
        :C => PkgData([:x1], NO_DEPS, NO_COMP),
    )
    info = pkg_info(data, keys(data); filter = false)
    Q = Problem([:A]; compat = Dict(:A => [:v3, :v2, :v1]))
    univ = prepare_pkg_info(info, Q)
    @test univ.info[:A].members == [[2, 3], [1, 4]] # the :C class laid out first
    @test univ.reps[:A] == [2, 4]
    @test resolve(info, Q) == Dict(:A => :v3, :C => :x1)

    R = Problem([:A])
    @test reranks(univ, R)
    @test resolve(info, R) == Dict(:A => :v4, :B => :w1)
    @test resolve_relaxed(univ, R) == Dict(:A => :v4, :B => :w1)
end

@testset "a relaxation names the version R admits, not the one Q did" begin
    # :A's :v3 and :v2 are indistinguishable to the registry — one class of two
    # members — and :v1 is a class of its own because it alone depends on :B.
    # Q's bound empties the first class outright; R lifts it, which both revives
    # the class and gives it a representative :A had none of under Q
    data = Dict(
        :A => PkgData([:v3, :v2, :v1],
            Dict(:v3 => Symbol[], :v2 => Symbol[], :v1 => [:B]), NO_COMP),
        :B => PkgData([:w1], NO_DEPS, NO_COMP),
    )
    info = pkg_info(data, keys(data); filter = false)
    Q = Problem([:A]; compat = Dict(:A => [:v1]))
    univ = prepare_pkg_info(info, Q)
    # the emptied class is still there, and still has no representative
    @test nclasses(univ.info[:A]) == 2
    @test univ.reps[:A] == [0, 3]
    @test resolve(info, Q) == Dict(:A => :v1, :B => :w1)

    R = Problem([:A])
    @test resolve(info, R) == Dict(:A => :v3)
    # naming through Q's representatives would index version 0 of :A; naming
    # through R's gives the version R's best class stands at
    @test resolve_relaxed(univ, R) == Dict(:A => :v3)

    # and with the bound merely loosened rather than dropped, the same class is
    # revived at a different member
    sat = SAT(univ)
    try
        @test resolve_relaxed(sat, Problem([:A])) == Dict(:A => :v3)
        @test resolve_relaxed(sat, Problem(Symbol[])) == Dict{Symbol,Symbol}()
        # Q itself is a relaxation of Q, and answers as Q answers
        @test resolve_relaxed(sat, Q) == resolve(info, Q)
    finally
        finalize(sat)
    end
end

@testset "push frames have to balance" begin
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]), NO_COMP),
        :B => PkgData([:w2, :w1], NO_DEPS, NO_COMP),
    )
    info = pkg_info(data, keys(data); filter = false)
    Q = Problem([:A]; compat = Dict(:A => [:v1]))
    sat = SAT(prepare_pkg_info(info, Q))
    try
        want = resolve_relaxed(sat, Problem([:A]))
        depth = sat.depth
        # asking from inside a temp-clause frame would have the swap pop that
        # frame instead of the query's, and nothing downstream could tell
        @test_throws AssertionError with_temp_clauses(sat) do
            resolve_relaxed(sat, Problem([:A]))
        end
        # the failed call left the instance alone, so it still answers
        @test sat.depth == depth
        @test resolve_relaxed(sat, Problem([:A])) == want
    finally
        finalize(sat)
    end
    # a body that pushes and does not pop is caught on the way out, by both
    # helpers, at the call that did it
    sat = SAT(prepare_pkg_info(info, Q))
    try
        @test_throws AssertionError with_temp_clauses(() -> sat_push(sat), sat)
    finally
        finalize(sat)
    end
    sat = SAT(prepare_pkg_info(info, Q))
    try
        @test_throws AssertionError with_classes_relaxed(() -> sat_push(sat), sat)
    finally
        finalize(sat)
    end
end

@testset "a ranking is a permutation stored as its deviation" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    # the identity is the empty deviation — the case every package a relaxation
    # left alone is in, and the one there is no branch for
    id = ClassRanking(collect(1:5))
    @test isempty(id.at) && isempty(id.to)
    @test collect(id) == 1:5
    @test collect(invperm(id)) == 1:5
    moved = 0
    for n = 1:6, _ = 1:20
        perm = randperm(n)
        r = ClassRanking(perm)
        # it is the permutation, read either way, through the `AbstractVector`
        # interface the descent takes a `view` of
        @test collect(r) == perm
        @test collect(invperm(r)) == invperm(perm)
        @test all(r[i] == perm[i] for i = 1:n)
        # and it stores only what moved
        @test length(r.at) == count(i -> perm[i] != i, 1:n)
        @test issorted(r.at)
        moved += length(r.at) > 0
    end
    @test moved > 0
end

@testset "one instance, many relaxations, no residue" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    for (m, n) in ((2, 3), (3, 2), (3, 3), (2, 4))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
        reqs = collect(make_reqs(rand(1:2^m-1)))
        info = pkg_info(data, keys(data); filter = false)
        compat, pins = random_constraints(m, n)
        bad = rand(1:n)
        Q = Problem(reqs; compat, pins,
            excludes = Pair{Symbol,Any}[:test => (p, v) -> v == bad])
        univ = prepare_pkg_info(info, Q)
        rels = relaxations(info, Q, 8)
        want = Dict(i => resolve(info, R) for (i, R) in enumerate(rels))
        sat = SAT(univ)
        try
            # `verdicts` (relaxation.jl) is the instance's answers to every
            # satisfiability question that can be put to it in terms of
            # packages and classes: two states that agree on it are
            # indistinguishable to anything above. The lifted reading is what
            # says the *frames* balance from the outside, `sat.depth` saying it
            # from the inside: `with_classes_relaxed` pops exactly one frame, so
            # one left open would have it lift that instead of the
            # deactivations, and the relaxed reading would come out wrong
            before = verdicts(sat)
            lifted = with_classes_relaxed(() -> verdicts(sat), sat)
            deact = copy(sat.deact)
            depth = sat.depth
            self = resolve(sat, Q.reqs)
            # each round asks the same relaxations in a fresh random order, so
            # a leak from one relaxation into the next would show up as an
            # answer that depends on what was asked before it
            for round = 1:4, i in randperm(length(rels))
                @test resolve_relaxed(sat, rels[i]) == want[i]
            end
            # the instance itself is as it was found: the same forbidden
            # classes, the same open frames, the same answers, and the same
            # answer to its own query — the last of which is the phase hints
            # too, since nothing but a solver knob distinguishes two descents
            # that add the same clauses
            @test sat.deact == deact
            @test sat.depth == depth
            @test verdicts(sat) == before
            @test with_classes_relaxed(() -> verdicts(sat), sat) == lifted
            @test resolve(sat, Q.reqs) == self
        finally
            finalize(sat)
        end
    end
end

@testset "the reachable set is downward closed in key_∅" begin
    # The walk keeps `{c : key_∅(c) ≤ t}` for a threshold that only rises. That
    # is what makes the kept set a prefix — in key_∅ order, not layout order —
    # and so representable by a pointer rather than a set. If it ever failed,
    # the walk would be keeping something its own rule does not license.
    checked = nontrivial = 0
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (2, 5))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:25
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            info = pkg_info(data, keys(data); filter = false)
            compat, pins = random_constraints(m, n)
            prob = Problem(reqs; compat, pins)
            univ = rank_pkg_info(info, prob)
            mark_installable!(univ.info)
            mark_reachable!(univ.info, prob.reqs, deactivations(univ), univ.ranks)
            for (p, ip) in univ.info
                k0 = univ.ranks[2][p]
                mc = nclasses(ip)
                kept = [j for j = 1:mc if ip.conflicts[j, end]]
                isempty(kept) && continue
                t = maximum(k0[j] for j in kept)
                @test all((k0[j] ≤ t) == ip.conflicts[j, end] for j = 1:mc)
                checked += 1
                length(kept) < mc && (nontrivial += 1)
            end
        end
    end
    # the sweep really did exercise packages that are not kept whole
    @test checked > 100
    @test nontrivial > 0
end
