"""
    version_permutations(info, order) :: Union{Nothing, Dict{P, Vector{Int}}}

Per package, the permutation that puts its versions into the ordering `order`
asks for: `perm[r]` is the index, in the T1 version list, of the `r`-th best
version. Only the packages that actually need reordering get an entry, and the
whole result is `nothing` when none of them do — so the *canonical* order, which
is the one the T1 artifact carries, costs nothing.

`order` is `nothing` for the canonical order, or a callable mapping a package to
a `lt` comparator over its versions (`lt(u, v)` = "u is preferred to v"). The
ordering is the one query knob that is not a *constraint*: it does not change
which solutions are valid, only which valid solution is optimal. That is why it
is a `resolve` parameter rather than part of the [`Problem`](@ref), and why the
T1 artifact can be built and cached without knowing it (see the manual's Theory
section, on order-free deletion and on interchangeability classes being
ordering-independent while the representative choice is not).
"""
version_permutations(
    info  :: AbstractDict{P, PkgInfo{P,V}},
    order :: Nothing,
) where {P,V} = nothing

function version_permutations(
    info  :: AbstractDict{P, PkgInfo{P,V}},
    order,
) where {P,V}
    perms = Dict{P, Vector{Int}}()
    for (p, info_p) in info
        vers = info_p.versions
        lt = order(p)
        # the common case by far — a comparator that agrees with the canonical
        # order, e.g. "newest first" on a registry that already lists versions
        # that way — is detected in one linear pass and costs nothing more
        issorted(vers; lt) && continue
        perms[p] = sortperm(vers; lt, alg = MergeSort)
    end
    return isempty(perms) ? nothing : perms
end

"""
    class_representatives(info, prob, perms) :: Dict{P, BitVector}

Per package, a mask of the versions that survive collapsing its
interchangeability classes: the T1 classes (`version_classes`) refined by
`prob`, then one representative per refined class.

**Refinement.** User compat and pins are the only constraints that can tell two
members of a T1 class apart, because they are stated on version *values*
instead of on the registry's rows — so splitting each class by `prob`'s
[`exclusion_masks`](@ref Resolver.exclusion_masks) is enough, and it costs
nothing for the packages the user did not constrain. Read as the virtual
package of `Problem`'s docstring, this is just row equality again, over one
more conflict column.

**Representative.** The best member in the active rank order — the
lowest-indexed one when the order is canonical, and otherwise the one `perms`
(see [`version_permutations`](@ref Resolver.version_permutations)) ranks first.

Collapsing is answer-preserving: identical constraint rows are in particular a
subset of each other, so a class member is dominated by its representative and
the redundancy theorem's swap argument applies (see the manual's Theory section,
lemma on interchangeable versions). Representatives *are* versions, and the
members dropped are exactly the ones the layered answer never chooses, so
nothing has to be mapped back afterwards.
"""
function class_representatives(
    info  :: AbstractDict{P, PkgInfo{P,V}},
    prob  :: Problem{P} = Problem(P[]),
    perms :: Union{Nothing, AbstractDict{P, Vector{Int}}} = nothing,
) where {P,V}
    excl = exclusion_masks(info, prob)
    keep = Dict{P,BitVector}()
    seen = Bool[] # refined classes already represented
    for (p, info_p) in info
        m = length(info_p.versions)
        cls = info_p.classes
        e = get(excl, p, nothing)
        perm = perms === nothing ? nothing : get(perms, p, nothing)
        k = falses(m)
        resize!(seen, 2m)
        fill!(seen, false)
        for r = 1:m
            # walk the versions best first, so the first member of a class
            # found is the one the active order ranks best
            i = perm === nothing ? r : perm[r]
            # refined class key: the T1 class, split by forbidden-ness
            c = 2 * cls[i] - (e !== nothing && e[i])
            seen[c] && continue
            seen[c] = true
            k[i] = true
        end
        keep[p] = k
    end
    return keep
end

"""
    prepare_pkg_info(info, prob, [info′]; group = true, order = nothing)

The per-resolve half of preprocessing, run on a T1 artifact (see
[`pkg_info`](@ref Resolver.pkg_info)) to produce the universe the SAT instance
is built over. Three steps, in this order:

1. **Lay out** each package's versions in the ordering the query wants
   (`version_permutations`), which is a no-op for the canonical order.
2. **Collapse** each interchangeability class, refined by `prob`, to its best
   member in that ordering (`class_representatives`).
3. **Filter** the collapsed universe for reachability and redundancy against
   the actual requirements and constraints (`filter_pkg_info!`).

Filtering has to come last: both of its passes are stated in terms of the
version ordering, and both get cheaper the fewer versions there are. The layout
and the collapse share one materialization with it: the collapse is handed over
as *marks* rather than as a rebuilt universe, so the filter's own first
`drop_unmarked!` performs both (see `copy_marked!`).

The result is a fresh dict — `info` is left untouched, so a T1 artifact stays
reusable across resolves — unless the caller passes its own scratch dict as
`info′`, or `info` itself when it owns it (which reordering rules out, since it
rebuilds the matrices). `group = false` skips the collapse, which is the
ungrouped path the tests compare against.
"""
function prepare_pkg_info(
    info  :: AbstractDict{P, PkgInfo{P,V}},
    prob  :: Problem{P},
    info′ :: Dict{P, PkgInfo{P,V}} = Dict{P, PkgInfo{P,V}}();
    group :: Bool = true,
    order = nothing, # nothing = canonical, else package -> `lt` comparator
) where {P,V}
    perms = version_permutations(info, order)
    keep = group ? class_representatives(info, prob, perms) : nothing
    # reordering reads every matrix while writing a differently laid out one,
    # so the in-place shortcut is only available for the canonical order
    perms === nothing || info′ !== info ||
        (info′ = Dict{P, PkgInfo{P,V}}())
    copy_marked!(info′, info, keep, perms)
    filter_pkg_info!(info′, prob)
    return info′
end

# Copy `info` into `info′`, laid out in the order `perms` gives (the order it
# already has when that is `nothing`) and marking active exactly the versions
# `keep` selects (all of them when it is `nothing`); when the two are the same
# dict, only the marks are written. This is how the collapse is handed to
# `filter_pkg_info!`: as *marks* rather than as an already-shrunken universe, so
# that the filter's own first `drop_unmarked!` materializes the collapse and its
# own deletions together, in one rebuild instead of two. Every pass the filter
# runs before that point reads the version flags as the version set, so they all
# see the collapsed problem — which is the problem the interchangeability lemma
# licenses substituting for the original.
function copy_marked!(
    info′ :: Dict{P, PkgInfo{P,V}},
    info  :: AbstractDict{P, PkgInfo{P,V}},
    keep  :: Union{Nothing, AbstractDict{P, BitVector}} = nothing,
    perms :: Union{Nothing, AbstractDict{P, Vector{Int}}} = nothing,
) where {P,V}
    src = Int[] # scratch: output column => source column
    for (p, info_p) in info
        X = info_p.conflicts
        m = length(info_p.versions)
        n = size(X, 2) - 1
        if info′ === info
            keep === nothing || (X[1:m, end] .= keep[p])
            continue
        end
        perm = perms === nothing ? nothing : get(perms, p, nothing)
        # a reordered partner's interaction block is indexed by that partner's
        # versions, so its columns move with them; every other column stays put
        cols = nothing
        if perms !== nothing
            for (q, off) in info_p.interacts
                perm_q = get(perms, q, nothing)
                perm_q === nothing && continue
                if cols === nothing
                    cols = resize!(src, n)
                    for j = 1:n
                        cols[j] = j
                    end
                end
                for (r, j) in enumerate(perm_q)
                    cols[off + r] = off + j
                end
            end
        end
        X′ = perm === nothing && cols === nothing ? copy(X) :
             relaid_conflicts(X, m, n, perm, cols)
        if keep === nothing
            X′[1:m, end] .= true
        elseif perm === nothing
            X′[1:m, end] .= keep[p]
        else
            k = keep[p]
            for r = 1:m
                X′[r, end] = k[perm[r]]
            end
        end
        X′[m+1, 1:n] .= true
        info′[p] = PkgInfo(
            perm === nothing ? copy(info_p.versions) : info_p.versions[perm],
            copy(info_p.depends), copy(info_p.interacts), X′,
            perm === nothing ? copy(info_p.classes) : info_p.classes[perm])
    end
    return info′
end

# A conflicts matrix laid out under a reordering: row `r` reads source row
# `perm[r]` (this package's versions) and column `j` reads source column
# `cols[j]` (a reordered partner's block); `nothing` for either means unchanged.
# Only the version rows of the `n` conflict columns are written — the caller
# sets the flag row and the version-flag column. Columns whose rows do not move
# are blitted word by word, which is the whole matrix whenever only partners
# were reordered.
function relaid_conflicts(
    X    :: BitMatrix,
    m    :: Int,
    n    :: Int,
    perm :: Union{Nothing, Vector{Int}},
    cols :: Union{Nothing, Vector{Int}},
)
    X′ = falses(size(X, 1), n + 1)
    W = col_words(X)
    ch = X.chunks
    ch′ = X′.chunks
    @inbounds for r = 1:n
        j = cols === nothing ? r : cols[r]
        base = (j - 1) * W
        base′ = (r - 1) * W
        if perm === nothing
            for w = 1:W
                ch′[base′ + w] = ch[base + w]
            end
        else
            for i′ = 1:m
                i = perm[i′]
                b = (ch[base + ((i - 1) >> 6) + 1] >> ((i - 1) & 63)) & 1
                ch′[base′ + ((i′ - 1) >> 6) + 1] |= b << ((i′ - 1) & 63)
            end
        end
    end
    return X′
end

filter_pkg_info!(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info),
) where {P,V} = filter_pkg_info!(info, Problem(reqs))

function filter_pkg_info!(
    info :: Dict{P, PkgInfo{P,V}},
    prob :: Problem{P},
) where {P,V}
    reqs = prob.reqs
    # user constraints enter as the exclusion masks of the virtual package
    # that represents them (see Problem.jl): reachability treats an excluded
    # version as one that conflicts with an always-present version, and
    # redundancy treats the mask as one more conflict column. the masks are
    # index-based, so they are rebuilt after every deletion round; only the
    # constrained packages cost anything
    excl = exclusion_masks(info, prob)
    # the version flags on entry are the version set to filter: normally all
    # of them, but `prepare_pkg_info` seeds them with the interchangeability
    # collapse, so that the first `drop_unmarked!` below deletes the collapsed
    # members along with everything the passes here find.
    #
    # arc consistency first: it needs no reachability marks and it is what
    # makes the deletions safe — dropping a package while a kept version
    # still depends on it would leave a dangling name. then redundancy
    # elimination — it needs no reachability marks and does most of the
    # shrinking — then alternate reachability and redundancy until neither
    # deletes anything: each deleted version can expose more of both kinds
    # of pruning, and every round preserves the resolved answer (see the
    # theory docs on iterating the filter). requirement packages survive
    # filtering unless they have no installable version at all — nothing
    # can satisfy those — so the requirements remain valid across rounds.
    # rounds strictly shrink the total version count, so the loop
    # terminates
    mark_installable!(info)
    mark_necessary!(info, excl)
    drop_unmarked!(info)
    while true
        total = sum(length(i.versions) for i in values(info); init = 0)
        excl = exclusion_masks(info, prob)
        mark_reachable!(info, reqs, excl)
        mark_necessary!(info, excl)
        drop_unmarked!(info)
        sum(length(i.versions) for i in values(info); init = 0) < total ||
            break
    end
    return info
end

"""
    find_reachable(info, reqs) :: Dict{P, Int}

This function finds a minimal "reachable" subset of package and versions that
could appear in pareto-optimal solutions to version resolution for the given set
of required "root" packages, using the following recursive logic:

- P in reqs => P[1] reachable
- P[i] reachable & P[i] depends on D => D[1] reachable
- P[i] reachable & P[i] conflicts w. reachable => P[i+1] reachable
- D[end] conflicts w. reachable & P[i] depends on D => P[i+1] reachable

The optional `excl` argument gives, per package, a mask of versions the user
forbade. Those are the conflict rows of the always-present virtual package that
represents the user constraints, so a reachable excluded version fires the same
degradation rule as any other conflict:

- P[i] reachable & P[i] excluded => P[i+1] reachable

The function returns a dictionary mapping packages to the maximum version index
of that package that could be reached in an optimal solution. If a package
cannot appear in an optimal solution, it will not appear in this dictionary.
"""
function find_reachable(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info),
    excl :: AbstractDict{P, BitVector} = EmptyDict{P,BitVector}(),
) where {P,V}
    # intern packages as integer indices so that the fixpoint loop below
    # hashes no package names; all derived tables are indexed by pkg id
    pkgs = sort!(collect(keys(info)))
    N = length(pkgs)
    ix = Dict{P,Int}(p => i for (i, p) in enumerate(pkgs))
    infos = Vector{PkgInfo{P,V}}(undef, N)
    nvers = Vector{Int}(undef, N)
    # interned depends, same order (id 0: dependency absent from info,
    # i.e. a package with no installable version — see below)
    deps = Vector{Vector{Int}}(undef, N)
    partners = Vector{Vector{Tuple{Int,Int}}}(undef, N)
    # interned exclusion masks (nothing: package unconstrained, the norm)
    excls = Vector{Union{Nothing,BitVector}}(nothing, N)
    for p = 1:N
        info_p = info[pkgs[p]]
        infos[p] = info_p
        nvers[p] = length(info_p.versions)
        excls[p] = get(excl, pkgs[p], nothing)
        deps[p] = Int[get(ix, q, 0) for q in info_p.depends]
        # (partner id, offset of p's version block in the partner's matrix)
        prt = Tuple{Int,Int}[
            (ix[q], info[q].interacts[pkgs[p]])
            for q in keys(info_p.interacts)
        ]
        partners[p] = sort!(prt)
    end

    # meaning (both map packages to version indices):
    #   - reach tracks fully processed reachable versions
    #   - queue tracks newly reachable versions not yet processed
    reach = zeros(Int, N)
    queue = zeros(Int, N)
    stack = Int[] # packages with queued work
    instack = falses(N)

    function enqueue(p::Int, j::Int)
        queue[p] = j
        if !instack[p]
            instack[p] = true
            push!(stack, p)
        end
    end

    # add next active version of p *after* i to the queue
    # do nothing if there's already version > i in reach/queue
    function next(p::Int, i::Int)
        reach[p] > i && return false
        queue[p] > i && return false
        m = nvers[p]
        X = infos[p].conflicts
        # first active version after i (the flag column is the active set)
        j = col_min_from(X, size(X, 2), i + 1, m)
        # j == 0: we're out of versions (i.e. saturated; see below)
        enqueue(p, j == 0 ? m + 1 : j)
        return true
    end

    rdeps = [Dict{Int,Int}() for _ = 1:N] # reverse dependency map
    # rdeps[p][q] == k means
    #   "k is latest reachable version of q that depends on p"

    for p in reqs
        # a requirement with no installable version reaches nothing
        i = get(ix, p, 0)
        i == 0 || enqueue(i, 1)
    end

    # notation:
    #   - p, q: packages
    #   - info_p = infos[p]
    #   - indices of versions of package p: i, j
    #   - indices of versions of package q: k
    while !isempty(stack)
        # get unprocessed package + version
        p = pop!(stack)
        instack[p] = false
        i = queue[p]
        queue[p] = 0
        info_p = infos[p]
        m = nvers[p]
        # check for saturation
        #   p saturated means: conflicts can force p to be uninstallable
        #   saturation represented by i > nvers[p]
        if i > m
            # p has become saturated
            for (q, k) in rdeps[p]
                # q@k depends on p, therefore
                # q@k conflicts with p being uninstallable
                # p being saturated means that can happen
                next(q, k)
            end
        end
        # process each newly reachable version of p
        excl_p = excls[p]
        for j = reach[p]+1:min(i, m)
            # user constraints: p@j conflicts with the always-present
            # version of the virtual package, so it can never be a resting
            # point — p can only be installed past it
            excl_p === nothing || !excl_p[j] || next(p, j)
            # dependencies
            for (k, q) in enumerate(deps[p])
                info_p.conflicts[j, k] || continue
                # p@j depends on q
                if q == 0
                    # q has no installable version at all, so neither has
                    # p@j: p can only be installed past it
                    next(p, j)
                    continue
                end
                rdeps_q = rdeps[q]
                rdeps_q[p] = max(get(rdeps_q, p, 0), j)
                next(q, 0) # q can be required
                # check if q is saturated:
                if reach[q] > nvers[q]
                    # p@j depends on q, therefore
                    # p@j conflicts with q being uninstallable
                    # q being saturated means that can happen
                    next(p, j)
                end
            end
            # find p@j's conflicts with reached versions of each partner:
            # conflicts are stored symmetrically, so scan the partner-side
            # column Y[k, c+j] == X[j, b+k], which is contiguous
            for (q, c) in partners[p]
                r = min(reach[q], nvers[q])
                r > 0 || continue
                k = col_max_upto(infos[q].conflicts, c + j, r)
                k == 0 && continue
                # p@j conflicts with q@k (the highest such k; pushing past
                # it subsumes pushing past any lower conflicting version)
                next(p, j)
                next(q, k)
            end
        end
        # update the reach map
        reach[p] = i
    end

    return Dict{P,Int}(pkgs[p] => reach[p] for p = 1:N if reach[p] > 0)
end

function mark_reachable!(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info),
    excl :: AbstractDict{P, BitVector} = EmptyDict{P,BitVector}(),
) where {P,V}
    reach = find_reachable(info, reqs, excl)
    for (p, info_p) in info
        r = min(get(reach, p, 0), length(info_p.versions))
        info_p.conflicts[1:r, end] .= true
        info_p.conflicts[r+1:end, end] .= false
    end
    return reach
end

function mark_installable!(
    info :: Dict{P, <: PkgInfo{P}},
) where {P}
    # a package with no active versions cannot be installed, so no version
    # depending on it can be either: deactivate those, to a fixpoint (each
    # deactivation can empty another package's active set). this is the
    # arc-consistency test — a sound approximation of deleting model-free
    # versions (see the theory docs) — and it is also what keeps `depends`
    # from naming a package that `drop_unmarked!` deletes as empty
    #
    # a dependency naming a package that is not in `info` at all counts the
    # same way: it cannot be satisfied either, and a deletion elsewhere (the
    # class collapse, say, or an earlier `drop_unmarked!`) is exactly how one
    # comes about. `info` does not change here, so the absentees are found once
    absent = Set{P}()
    for info_p in values(info), q in info_p.depends
        haskey(info, q) || push!(absent, q)
    end
    empties = Set{P}()
    while true
        empty!(empties)
        union!(empties, absent)
        for (p, info_p) in info
            X = info_p.conflicts
            any(X[i, end] for i = 1:length(info_p.versions)) && continue
            push!(empties, p)
        end
        isempty(empties) && return info
        dirty = false
        for info_p in values(info)
            X = info_p.conflicts
            for (k, q) in enumerate(info_p.depends)
                q in empties || continue
                # q is uninstallable, drop the p@i that depend on it
                for i = 1:length(info_p.versions)
                    X[i, end] || continue
                    X[i, k] || continue
                    X[i, end] = false
                    dirty = true
                end
            end
        end
        dirty || return info
    end
end

# one sweep of the domination candidates against an exclusion column: every
# still-tracked version the column excludes has its candidates restricted to
# the versions the column also excludes, and drops out of the sweep when none
# are left. A named function rather than a closure in `mark_necessary!`'s
# per-package loop, which is the filter's hot path.
function exclude_column!(
    D    :: Vector{UInt64}, # per-version candidate masks
    T    :: Vector{UInt64}, # versions still tracked
    E    :: Vector{UInt64}, # scratch: the column, padded to the column width
    mask :: BitVector,
    W    :: Int,
)
    resize!(E, W)
    fill!(E, 0)
    copyto!(E, 1, mask.chunks, 1, min(W, length(mask.chunks)))
    @inbounds for w = 1:W
        c = E[w] & T[w]
        while !iszero(c)
            i = ((w - 1) << 6) + trailing_zeros(c) + 1
            c &= c - 1
            o = (i - 1) * W
            live = UInt64(0)
            for w′ = 1:W
                live |= (D[o + w′] &= E[w′])
            end
            iszero(live) && (T[w] &= ~(UInt64(1) << ((i - 1) & 63)))
        end
    end
    return T
end

function mark_necessary!(
    info :: Dict{P, PkgInfo{P,V}},
    excl :: AbstractDict{P, BitVector} = EmptyDict{P,BitVector}(),
) where {P,V}
    # intern packages as integer indices so that the work loop below
    # hashes no package names (sorted so the processing order — and with
    # it any order-dependent tie — is deterministic)
    pkgs = sort!(collect(keys(info)))
    N = length(pkgs)
    ix = Dict{P,Int}(p => i for (i, p) in enumerate(pkgs))
    infos = Vector{PkgInfo{P,V}}(undef, N)
    nvers = Vector{Int}(undef, N)
    partners = Vector{Vector{NTuple{3,Int}}}(undef, N)
    # interned exclusion masks (nothing: package unconstrained, the norm)
    excls = Vector{Union{Nothing,BitVector}}(nothing, N)
    for p = 1:N
        info_p = info[pkgs[p]]
        infos[p] = info_p
        nvers[p] = length(info_p.versions)
        excls[p] = get(excl, pkgs[p], nothing)
        # (partner id, partner's version block offset in p's matrix,
        #  p's version block offset in the partner's matrix)
        prt = NTuple{3,Int}[
            (ix[q], b, info[q].interacts[pkgs[p]])
            for (q, b) in info_p.interacts
        ]
        partners[p] = sort!(prt)
    end
    # some work buffers
    A = UInt64[]        # active version mask
    D = UInt64[]        # per-version domination candidate masks
    T = UInt64[]        # mask of versions still tracked by the sweep
    E = UInt64[]        # exclusion mask, padded to the column width
    R = Int[]           # redundant indices vector
    # initialize active column flags
    for p = 1:N
        info_p = infos[p]
        X = info_p.conflicts
        m = nvers[p]
        n = size(X, 2) - 1
        resize!(A, col_words(X))
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        for j = 1:n
            # we don't have to look at columns that have no
            # conflicts for any active versions of package
            X[m+1, j] = col_intersects(X, j, A)
        end
        # conflict columns for inactive partner versions are also
        # irrelevant: no model over the active versions can include
        # the partner, so the conflict cannot constrain anything
        for (q, b, c) in partners[p]
            Y = infos[q].conflicts
            for k = 1:nvers[q]
                X[m+1, b+k] &= Y[k, end]
            end
        end
    end
    # process every package once, cheapest first (stable sort breaks ties
    # by name); packages woken by a partner's deletions are appended to
    # the queue unless already pending
    queue = sort!(collect(1:N), by = p -> length(partners[p]))
    inwork = trues(N)
    head = 1
    while head ≤ length(queue)
        p = queue[head]
        head += 1
        inwork[p] || continue
        inwork[p] = false
        # get conflicts & dimensions
        info_p = infos[p]
        X = info_p.conflicts
        m = nvers[p]
        m > 1 || continue # unique version cannot be redundant
        n = size(X, 2) - 1
        W = col_words(X)
        # active version mask
        resize!(A, W)
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        nact = 0
        @inbounds for w = 1:W
            nact += count_ones(A[w])
        end
        nact > 1 || continue # unique version cannot be redundant
        # find redundant versions: i < j and X[i, k] => X[j, k] for all
        # active k means an earlier version is strictly more compatible,
        # therefore i will always be chosen instead of j. the candidates
        # dominated by i are the active versions worse than i that have a
        # constraint in every column i does. sweep the active columns
        # once, restricting the candidates of every version with a bit in
        # the column to the column's row set; versions whose candidate
        # set empties drop out of the sweep. a version dominated only by
        # dominated versions is also dominated by their (transitively
        # live) dominators, so batching over the initial active set finds
        # exactly the sequentially dominated versions.
        resize!(D, m * W)
        resize!(T, W)
        copyto!(T, A)
        @inbounds for w = 1:W
            c = A[w]
            while !iszero(c)
                i = ((w - 1) << 6) + trailing_zeros(c) + 1
                c &= c - 1
                # candidates of i: active versions worse than i
                o = (i - 1) * W
                for w′ = 1:W
                    D[o + w′] = A[w′]
                end
                wt = i >> 6
                for w′ = 1:min(wt, W)
                    D[o + w′] = 0
                end
                r = i & 63
                if r != 0 && wt < W
                    D[o + wt + 1] &= ~((UInt64(1) << r) - 1)
                end
                live = UInt64(0)
                for w′ = 1:W
                    live |= D[o + w′]
                end
                iszero(live) && (T[w] &= ~(UInt64(1) << ((i - 1) & 63)))
            end
        end
        prev = -1 # chunk base of the previous active column
        for k = 1:n
            X[m+1, k] || continue # inactive column
            base = (k - 1) * W
            # compat bounds are ranges, so identical adjacent columns are
            # common — and re-ANDing an identical column changes nothing
            if prev ≥ 0
                same = true
                @inbounds for w = 1:W
                    if X.chunks[base + w] != X.chunks[prev + w]
                        same = false
                        break
                    end
                end
                same && continue
            end
            prev = base
            @inbounds for w = 1:W
                c = X.chunks[base + w] & T[w]
                while !iszero(c)
                    i = ((w - 1) << 6) + trailing_zeros(c) + 1
                    c &= c - 1
                    o = (i - 1) * W
                    live = UInt64(0)
                    for w′ = 1:W
                        live |= (D[o + w′] &= X.chunks[base + w′])
                    end
                    iszero(live) && (T[w] &= ~(UInt64(1) << ((i - 1) & 63)))
                end
            end
        end
        # the virtual package representing the constraints contributes one
        # more conflict column -- the exclusion mask. its version is always
        # present, so the column is always active: an excluded version must not
        # dominate a non-excluded one
        excl_p = excls[p]
        excl_p === nothing || exclude_column!(D, T, E, excl_p, W)
        # union of all candidate sets = the dominated versions
        fill!(T, UInt64(0)) # reuse as the union accumulator
        @inbounds for w = 1:W
            c = A[w]
            while !iszero(c)
                i = ((w - 1) << 6) + trailing_zeros(c) + 1
                c &= c - 1
                o = (i - 1) * W
                for w′ = 1:W
                    T[w′] |= D[o + w′]
                end
            end
        end
        empty!(R)
        @inbounds for w = 1:W
            c = T[w]
            while !iszero(c)
                push!(R, ((w - 1) << 6) + trailing_zeros(c) + 1)
                c &= c - 1
            end
        end
        isempty(R) && continue
        # deactivate redundant versions
        X[R, end] .= false
        for (q, b, c) in partners[p]
            infos[q].conflicts[nvers[q] + 1, c .+ R] .= false
            if !inwork[q] # can create new redundancies
                inwork[q] = true
                push!(queue, q)
            end
        end
    end
end

function drop_unmarked!(
    info′ :: Dict{P, <: PkgInfo{P}},
    info  :: Dict{P, <: PkgInfo{P}} = info′,
) where {P}
    # first pass: per package, compute the kept versions and kept columns
    # from the (definitive) version flags — all of it before any package
    # is rebuilt, since kept columns are read off the partners' matrices
    masks = Dict{P, Tuple{BitVector, BitVector}}()
    A = UInt64[] # active version mask buffer
    for (p, info_p) in info
        X = info_p.conflicts
        m = length(info_p.versions)
        n = size(X, 2) - 1
        W = col_words(X)
        # active versions
        I = X[1:m, end]
        resize!(A, W)
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        # kept columns: dependency columns some active version depends
        # on, and conflict columns of active partner versions — a
        # partner's column block must stay aligned with its version list,
        # so kept partner versions keep their columns even if empty
        K = falses(n)
        for k = 1:length(info_p.depends)
            K[k] = col_intersects(X, k, A)
        end
        for (q, b) in info_p.interacts
            Y = info[q].conflicts
            copy_col_bits!(K.chunks, b, Y, size(Y, 2), length(info[q].versions))
        end
        masks[p] = (I, K)
    end
    # save original version counts (needed if info′ === info)
    N = Dict{P,Int}(
        p => length(info_p.versions)
        for (p, info_p) in info
    )
    # second pass: shrink each PkgInfo
    for (p, info_p) in info
        # abbreviate components
        V = info_p.versions
        D = info_p.depends
        T = info_p.interacts
        X = info_p.conflicts
        C = info_p.classes
        m = length(V)
        n = size(X, 2) - 1
        W = col_words(X)
        I, K = masks[p]
        m′ = count(I)
        # delete if no active versions
        if m′ == 0
            delete!(info′, p)
            continue
        end
        # keep as is if everything is active (with the column flags
        # normalized, as rebuilding would leave them)
        if m′ == m && all(K)
            X[m+1, 1:n] .= true
            if info !== info′
                info′[p] = PkgInfo(copy(V), copy(D), copy(T), copy(X), copy(C))
            end
            continue
        end
        # compute shrunken components
        V′ = V[I]
        D′ = D[K[1:length(D)]]
        T′ = Dict{P,Int}()
        b′ = length(D′)
        for (q, b) in sort!(collect(T), by=last)
            n′ = count(K[b .+ (1:N[q])])
            if n′ > 0
                T′[q] = b′
                b′ += n′
            end
        end
        R′ = padded_rows(m′)
        W′ = R′ >> 6
        X′ = falses(R′, b′ + 1)
        src = X.chunks
        dst = X′.chunks
        # kept-version mask words for the gather below
        resize!(A, W)
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        prefix = findlast(I) == m′ # kept versions are 1:m′
        j′ = 0
        for j = 1:n
            K[j] || continue
            j′ += 1
            sb = (j - 1) * W
            db = (j′ - 1) * W′
            if prefix
                # kept versions are a prefix: straight word copy
                nw = (m′ + 63) >> 6
                @inbounds for w = 1:nw
                    dst[db + w] = src[sb + w]
                end
                r = m′ & 63
                r != 0 && @inbounds (dst[db + nw] &= (UInt64(1) << r) - 1)
            else
                # gather the kept versions' bits
                acc = UInt64(0)
                na = 0
                dw = db + 1
                @inbounds for w = 1:W
                    avail = A[w]
                    iszero(avail) && continue
                    v = src[sb + w]
                    while !iszero(avail)
                        acc |= ((v >> trailing_zeros(avail)) & 1) << na
                        avail &= avail - 1
                        na += 1
                        if na == 64
                            dst[dw] = acc
                            dw += 1
                            acc = UInt64(0)
                            na = 0
                        end
                    end
                end
                na > 0 && @inbounds (dst[dw] = acc)
            end
        end
        @assert j′ == b′
        X′[1:m′, end] .= true  # kept versions are active
        X′[m′+1, 1:b′] .= true # kept columns are active
        # the surviving versions keep the classes they were in: dropping
        # versions and columns can only *merge* row-equality classes (fewer
        # rows to distinguish, fewer columns to differ in), so restricting the
        # partition stays sound, if possibly finer than the truth. `pkg_info`
        # is where it is computed exactly
        C′ = renumber_classes!(C[I])
        # assign new struct into info′
        info′[p] = PkgInfo(V′, D′, T′, X′, C′)
    end
    return info′
end

function drop_excluded!(
    info :: Dict{P, PkgInfo{P,V}},
    drop :: SetOrVec{P},
) where {P,V}
    drop isa AbstractSet || (drop = Set{P}(drop))
    for (p, info_p) in info
        # deactivate all versions of dropped packages
        p ∈ drop || continue
        info_p.conflicts[:, end] .= false
    end
    # deactivate versions that depend on dropped packages
    mark_installable!(info)
    drop_unmarked!(info)
end

function check_info_structure(
    info :: Dict{P, PkgInfo{P,V}},
) where {P,V}
    for (p, info_p) in info
        for q in info_p.depends
            q in keys(info) ||
                return :depends, p, q
        end
        interacts = sort!(collect(info_p.interacts), by=last)
        for (i, (q, b)) in enumerate(interacts)
            if i == 1
                b == length(info_p.depends) ||
                    return :conflicts, p, q
            end
            if i < length(interacts)
                b + length(info[q].versions) == interacts[i+1][2] ||
                    return :conflicts, p, q
            end
        end
    end
end
