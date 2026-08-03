# Phase B of a diagnosis: bound-level detail, on per-cluster closure
# sub-instances.
#
# The production instance keeps third-party bounds as hard clauses -- selector
# overhead on the hot path is not worth paying for a case that almost never
# happens -- so relaxing one needs a separate instance. Two results make that
# cheap and exact:
#
#   * Lemma D3 (closure exactness): a query assuming only a cluster's
#     requirements is decided by the sub-problem on their dependency closure,
#     so the sub-instance can be tiny compared to the universe.
#
#   * Theorem D1: the closure of the *filtered* universe answers as the raw
#     problem would -- but only for column-closed groups. For conflicts that
#     means an unordered *pair* of packages, never a single compat
#     declaration: Proposition D1' is what happens at finer granularity. So a
#     bound-level story blames incompatible pairs, and an upstream suggestion
#     names both sides.
#
# The queries themselves are plain satisfiability checks on a rebuilt
# instance: with the closure restricted, an instance costs microseconds, and a
# cluster needs one per pair. That is far simpler than guarding every conflict
# with a selector, and the sub-instances are small by construction.

# Phase B rebuilds its sub-instance once per probe, and needs one probe per
# incompatibility to shrink the story plus one per survivor to test it against
# an upstream release. A probe costs about one SAT build over the closure, so
# the work goes as (incompatibilities × versions in the closure) -- measured at
# roughly 3 µs a unit on the General registry, where every closure carries
# `julia` and its couple of thousand versions along with it.
#
# The budget below is therefore about a second of probing, which is a fair
# price on a path that has already failed. Beyond it the report stops at
# requirement level rather than spend minutes: a conflict whose closure is that
# big belongs to a hub package, and those are nearly always explained by the
# user's own constraints anyway.
#
# (The ceiling is an artifact of rebuilding the instance per probe. Guarding
# each pair's conflicts with a selector, the way the production instance guards
# the user constraints, would make a probe one solve instead of one build and
# retire the budget entirely; it needs SAT.jl's rectangle encoder factored out
# to share, which is why it is not done here.)
const BOUND_STORY_BUDGET = 400_000

# every package reachable from `reqs` along *any* version's dependency edges
function dep_closure(info::Dict{P,<:PkgInfo{P}}, reqs) where {P}
    seen = Set{P}(p for p in reqs if haskey(info, p))
    queue = collect(seen)
    while !isempty(queue)
        p = pop!(queue)
        for q in info[p].depends
            q ∈ seen && continue
            push!(seen, q)
            push!(queue, q)
        end
    end
    return seen
end

# `info` restricted to `keep`: same versions and dependencies (the closure is
# dependency-closed, so none dangle), conflicts with packages outside dropped.
# Lemma D3 says that loses nothing -- an outside package is never forced, and
# a model of the restriction extends by leaving it out.
function restrict_info(
    info :: Dict{P,PkgInfo{P,V}},
    keep :: AbstractSet{P},
) where {P,V}
    out = Dict{P,PkgInfo{P,V}}()
    for p in keep
        info_p = info[p]
        m = length(info_p.versions)
        D = copy(info_p.depends)
        parts = sort!(P[q for q in keys(info_p.interacts) if q ∈ keep])
        T = Dict{P,Int}()
        off = length(D)
        for q in parts
            T[q] = off
            off += length(info[q].versions)
        end
        X = falses(padded_rows(m), off + 1)
        Y = info_p.conflicts
        for k = 1:length(D), i = 1:m
            X[i, k] = Y[i, k]
        end
        for q in parts
            b = info_p.interacts[q]
            c = T[q]
            for j = 1:length(info[q].versions), i = 1:m
                X[i, c+j] = Y[i, b+j]
            end
        end
        X[1:m, end] .= true   # every version active
        X[m+1, 1:off] .= true # every column active
        # every version survives and only columns go, and dropping columns can
        # only *merge* row-equality classes, so carrying the partition over is
        # sound (possibly finer than the truth, which is all `classes` promises)
        out[p] = PkgInfo{P,V}(copy(info_p.versions), D, T, X,
                              copy(info_p.classes))
    end
    return out
end

# the unordered pairs of packages with a conflict between them: the finest
# column-closed conflict groups there are
function interacting_pairs(info::Dict{P,<:PkgInfo{P}}) where {P}
    pairs = Tuple{P,P}[]
    for (p, info_p) in info, q in keys(info_p.interacts)
        p < q || continue
        m = length(info_p.versions)
        b = info_p.interacts[q]
        n = length(info[q].versions)
        any(info_p.conflicts[i, b+j] for i = 1:m, j = 1:n) || continue
        push!(pairs, (p, q))
    end
    return sort!(pairs)
end

# restore (`on`) or clear (`off`) pair (p, q)'s conflict block in `work`,
# reading the original bits from the pristine `src`; the two have the same
# layout, `work` being a copy
function set_pair!(
    work :: Dict{P,PkgInfo{P,V}},
    src  :: Dict{P,PkgInfo{P,V}},
    p    :: P,
    q    :: P,
    on   :: Bool,
) where {P,V}
    for (a, b) in ((p, q), (q, p))
        Xw = work[a].conflicts
        Xs = src[a].conflicts
        off = work[a].interacts[b]
        for j = 1:length(work[b].versions), i = 1:length(work[a].versions)
            Xw[i, off+j] = on && Xs[i, off+j]
        end
    end
end

copy_info(info::Dict{P,PkgInfo{P,V}}) where {P,V} =
    Dict{P,PkgInfo{P,V}}(
        p => PkgInfo{P,V}(copy(i.versions), copy(i.depends),
                          copy(i.interacts), copy(i.conflicts), copy(i.classes))
        for (p, i) in info)

# `prob` as it applies inside a closure sub-instance, with the constraints of
# the packages in `on` in force and everything else's relaxed. Admission knobs
# have to be masked rather than dropped: a knob is stated about versions, not
# about a package, so switching it off for one package means excusing that
# package from its predicate.
function closure_problem(
    prob :: Problem{P},
    reqs :: AbstractVector{P},
    on   :: AbstractSet{P},
) where {P}
    is_constrained(prob) && !isempty(on) || return Problem(reqs)
    compat = Dict{P,valtype(prob.compat)}(
        p => s for (p, s) in prob.compat if p ∈ on)
    pins = Dict{P,valtype(prob.pins)}(
        p => v for (p, v) in prob.pins if p ∈ on)
    excludes = Pair{Symbol,Any}[
        kind => ((p, v) -> p ∈ on && forbids(p, v)::Bool)
        for (kind, forbids) in prob.excludes]
    Problem(reqs; compat, pins, excludes)
end

# A pair's conflicts, as facts: group one side's versions by which versions of
# the other they rule out, so the report reads as a handful of statements rather
# than a cell-by-cell dump.
#
# Most of those statements are about versions the query already excludes. A
# minimal conflict needs them -- something has to close off the old versions of
# a package, or the conflict would not be a conflict -- but a reader wants the
# one line about versions they could actually get, not five siblings at the same
# level. So each group is marked `incidental` when every version of `p` in it is
# already forbidden, and reports demote those. If *every* group is incidental
# nothing is demoted: better a verbose explanation than none.
function pair_facts(
    info :: Dict{P,PkgInfo{P,V}},
    prob :: Problem{P},
    p    :: P,
    q    :: P,
) where {P,V}
    info_p = info[p]
    vers_q = info[q].versions
    b = info_p.interacts[q]
    m = length(info_p.versions)
    n = length(vers_q)
    groups = Dict{Vector{Int},Vector{Int}}() # pattern over q => versions of p
    for i = 1:m
        pat = Int[j for j = 1:n if info_p.conflicts[i, b+j]]
        isempty(pat) && continue
        push!(get!(() -> Int[], groups, pat), i)
    end
    pats = sort!(collect(keys(groups)))
    live = Dict{Vector{Int},Bool}(
        pat => any(!is_excluded(prob, p, info_p.versions[i])
                   for i in groups[pat])
        for pat in pats)
    demote = any(values(live))
    facts = Fact[]
    for pat in pats
        bad = Set(pat)
        push!(facts, Bound{P,V}(
            p, V[info_p.versions[i] for i in groups[pat]],
            q, V[vers_q[j] for j = 1:n if j ∉ bad],
            demote && !live[pat]))
    end
    return facts
end

# the pair-level bound: which versions of `p` conflict with something of `q`
# at all, and which versions of `q` survive all of them
function pair_bound(info::Dict{P,PkgInfo{P,V}}, p::P, q::P) where {P,V}
    info_p = info[p]
    vers_q = info[q].versions
    b = info_p.interacts[q]
    m = length(info_p.versions)
    n = length(vers_q)
    pv = V[info_p.versions[i] for i = 1:m
           if any(info_p.conflicts[i, b+j] for j = 1:n)]
    qa = V[vers_q[j] for j = 1:n
           if !any(info_p.conflicts[i, b+j] for i = 1:m)]
    Bound{P,V}(p, pv, q, qa)
end

"""
    bound_story(info, prob, creqs, users, by) -> (facts, upstream, kept)

Bound-level detail for one conflict cluster: which package-pair
incompatibilities and which of the user's own constraints are needed to make
`creqs` unsatisfiable, and which incompatibilities a single upstream release
could dissolve.

Runs on a sub-instance restricted to `creqs`' dependency closure, relaxing
conflicts one *pair* at a time — the finest granularity Theorem D1 covers. The
deletion order is actionability-biased: incompatibilities only an upstream
release can move go first, so a conflict tellable either way keeps the
explanation the user can act on. `kept` is the subset of `users` that survived;
it is `nothing` when the closure was too big to be worth the solves, and the
caller has to shrink the user groups itself.
"""
function bound_story(
    info  :: Dict{P,PkgInfo{P,V}},
    prob  :: Problem{P},
    creqs :: Vector{P},
    users :: Vector{P},
    by    :: Function,
) where {P,V}
    none = (Fact[], UpstreamFix{P,V}[], nothing)
    keep = dep_closure(info, creqs)
    isempty(keep) && return none
    src = restrict_info(info, keep)
    userset = Set{P}(p for p in users if p ∈ keep)
    # filter the sub-problem for *this cluster's* requirements before probing
    # it. The universe was prepared for all of them at once, and a cluster is
    # usually a couple of packages, so reachability prunes hard here -- which
    # matters a lot, since every probe below rebuilds an instance over whatever
    # is left. Filtering with a subset of the requirements is licensed by the
    # requirement-monotonicity proposition, and doing it with every group in
    # force is what puts the result back under Theorem D1.
    filter_pkg_info!(src, closure_problem(prob, creqs, userset))
    all(p -> haskey(src, p), creqs) || return none
    pairs = interacting_pairs(src)
    nvers = sum(length(i.versions) for i in values(src); init = 0)
    length(pairs) * nvers ≤ BOUND_STORY_BUDGET || return none
    work = copy_info(src)

    # is the cluster satisfiable with these pairs in force? leaves `work`
    # holding exactly the relaxation it tested, so the caller can resolve on it.
    # `subprob` is hoisted because it is constant across a deletion pass and
    # rebuilding its constraint dictionaries per probe is pure overhead
    function probe(on::AbstractVector{Bool}, subprob::Problem{P})
        for (k, (p, q)) in enumerate(pairs)
            set_pair!(work, src, p, q, on[k])
        end
        sat = SAT(work, subprob)
        try is_satisfiable(sat, creqs)
        finally
            finalize(sat)
        end
    end

    # Every constraint the query really has stays in force for the probes and
    # for the solutions they report: an upstream release is only worth
    # suggesting if it fixes the conflict under the constraints the user
    # actually gave, and a solution that quietly ignored one -- a julia bound,
    # say, or the ban on prereleases -- would be no solution at all. Shrinking
    # the user set is a separate question, about which facts the *story* needs.
    subprob = closure_problem(prob, creqs, userset)

    # everything on must reproduce the failure; closure exactness says it does
    on = trues(length(pairs))
    probe(on, subprob) && return none

    # group MUS by a single biased deletion pass: the incompatibilities first,
    # so a conflict the user's own constraints can explain is blamed on them
    for k in eachindex(pairs)
        on[k] = false
        probe(on, subprob) && (on[k] = true)
    end
    story = copy(userset)
    for p in sort!(collect(userset))
        delete!(story, p)
        probe(on, closure_problem(prob, creqs, story)) &&
            push!(story, p) # needed: keep it
    end

    facts = Fact[]
    for (k, (p, q)) in enumerate(pairs)
        on[k] || continue
        append!(facts, pair_facts(src, prob, p, q))
    end

    # single-pair probes: with every other incompatibility and every surviving
    # user constraint still in force, would dissolving just this one resolve
    # the cluster? if so, one upstream release on either side would do it.
    #
    ups = UpstreamFix{P,V}[]
    for (k, (p, q)) in enumerate(pairs)
        on[k] || continue
        on[k] = false
        if probe(on, subprob)
            sol = resolve_prepared(work, subprob; by, diagnose = false)
            sol === nothing ||
                push!(ups, UpstreamFix{P,V}(pair_bound(src, p, q), sol))
        end
        on[k] = true
    end
    kept = P[p for p in users if p ∈ story]
    return facts, rank_upstream!(ups, src), kept
end

"""
    rank_upstream!(ups, info) -> ups

Order upstream suggestions by how plausible they are to act on, most first.

Two keys, in order:

1. **How current the blamed versions are.** A release that relaxes a bound its
   *newest* versions carry is a thing a maintainer might actually do; asking
   them to reach back and loosen a range that only ancient versions had is not.
   Measured as the best rank the incompatibility touches on either side —
   versions are listed best-first, so smaller is newer.
2. **How good the result is.** Between two plausible asks, the one that gets you
   newer packages is the better suggestion. Measured on the witness, which is an
   optimal layered solution, over the two packages named.

Ties break on the package names, so the order is deterministic.
"""
function rank_upstream!(
    ups  :: Vector{UpstreamFix{P,V}},
    info :: Dict{P,PkgInfo{P,V}},
) where {P,V}
    rank(p, v) = something(findfirst(==(v), info[p].versions), typemax(Int))
    best(p, vs) = isempty(vs) ? typemax(Int) : minimum(rank(p, v) for v in vs)
    function key(u::UpstreamFix{P,V})
        b = u.bound
        blamed = min(best(b.pkg, b.versions),
                     best(b.dep, V[v for v in info[b.dep].versions
                                   if v ∉ b.allowed]))
        got = sum(haskey(u.solution, p) ? rank(p, u.solution[p]) : 0
                  for p in (b.pkg, b.dep))
        (blamed, got, b.pkg, b.dep)
    end
    return sort!(ups; by = key)
end
