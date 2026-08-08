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
# them. Two sharper statements are checked alongside, because a failure in
# either localizes the cause: that the Q-filtered universe contains everything
# a fresh filter for R would keep (Theorem B), and that the reachability walk's
# kept set is downward closed in `key_∅` (the Remark), which is the invariant
# that lets it carry an integer frontier instead of a set.

using Resolver: PkgData, PkgInfo, Problem, Universe, NoKinds, pkg_info, nclasses,
    prepare_pkg_info, rank_pkg_info, resolve, resolve_prepared, class_ranking,
    copy_ranked!, version_permutations, deactivations, mark_installable!,
    mark_reachable!

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

function relax(Q::Problem{P}, drop_reqs, drop_srcs) where {P}
    reqs = P[p for p in Q.reqs if p ∉ drop_reqs]
    nocomp = Set(p for (k, p) in drop_srcs if k === :compat)
    nopin  = Set(p for (k, p) in drop_srcs if k === :pin)
    nokind = Set(k for (t, k) in drop_srcs if t === :kind)
    compat = filter(kv -> first(kv) ∉ nocomp, Dict(Q.compat))
    pins   = filter(kv -> first(kv) ∉ nopin,  Dict(Q.pins))
    excludes = Pair{Symbol,Any}[kv for kv in Q.excludes if first(kv) ∉ nokind]
    return Problem(reqs; compat, pins, excludes)
end

# a sample: the bottom of the lattice (everything dropped) first, since it is
# the hardest case and the one caching depends on, then random points
function relaxations(Q::Problem{P}, k::Int) where {P}
    srcs = sources(Q)
    out = [relax(Q, P[], srcs)]
    for _ = 1:k
        dr = P[p for p in unique(Q.reqs) if rand() < 0.35]
        ds = Any[s for s in srcs if rand() < 0.5]
        isempty(dr) && isempty(ds) && continue
        push!(out, relax(Q, dr, ds))
    end
    return out
end

# What reuse means: keep the filtered universe, re-rank its classes for `R`
# — new representatives and, where `R` moves them, a new layout — and solve.
# No pass of the filter is run again.
function resolve_reuse(univ::Universe{P,V}, R::Problem{P}) where {P,V}
    info = univ.info
    reps, cperms = class_ranking(info, R, version_permutations(info, nothing))
    u = cperms === nothing ?
        Universe{P,V}(info, reps, nothing) :
        copy_ranked!(
            Universe{P,V}(Dict{P,PkgInfo{P,V}}(), Dict{P,Vector{Int}}(), nothing),
            info, reps, cperms)
    return resolve_prepared(u, R)
end

# the (package, version-tuple) identity of every class left in a universe, so
# that two independently built universes can be compared by content
classkeys(univ) = Set(
    (p, Tuple(ip.versions[j] for j in mem))
    for (p, ip) in univ.info for mem in ip.members)

@testset "a Q-filtered universe answers every relaxation of Q" begin
    reused = 0
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (2, 5))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:30
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            info = pkg_info(data, keys(data); filter = false)
            compat, pins = random_constraints(m, n)
            bad = rand(1:n)
            excludes = rand(Bool) ?
                Pair{Symbol,Any}[:test => (p, v) -> v == bad] : NoKinds
            Q = Problem(reqs; compat, pins, excludes)
            univ = prepare_pkg_info(info, Q)
            have = classkeys(univ)
            for R in relaxations(Q, 4)
                # Theorem C: reuse answers R exactly as a fresh resolve does
                @test resolve_reuse(univ, R) == resolve(info, R)
                # Theorem B, the sharper statement: nothing a fresh filter for
                # R keeps is missing from the universe filtered for Q
                @test isempty(setdiff(classkeys(prepare_pkg_info(info, R)), have))
                reused += 1
            end
        end
    end
    @test reused > 900
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
