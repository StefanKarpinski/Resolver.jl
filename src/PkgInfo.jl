struct PkgInfo{P,V}
    versions  :: Vector{V}
    depends   :: Vector{P}
    interacts :: Dict{P, Int}
    conflicts :: BitMatrix
    # the interchangeability partition of `versions`, as a class id per
    # version (see `version_classes`): versions share a class when nothing in
    # the registry can tell them apart. Derived from `conflicts`, and at the
    # T1 boundary (`pkg_info`) exactly so
    classes   :: Vector{Int}
end

PkgInfo(
    versions  :: Vector{V},
    depends   :: Vector{P},
    interacts :: Dict{P,Int},
    conflicts :: BitMatrix,
) where {P,V} = PkgInfo{P,V}(versions, depends, interacts, conflicts,
    version_classes(conflicts, length(versions)))

# `classes` is a function of `conflicts`, so it is not compared: in-place
# shrinking (see `drop_unmarked!`) leaves a sound but possibly finer partition
# than a recomputation would give, and that difference is not a difference of
# content
function Base.:(==)(a::PkgInfo, b::PkgInfo)
    a.versions  == b.versions  &&
    a.depends   == b.depends   &&
    a.interacts == b.interacts &&
    a.conflicts == b.conflicts
end

"""
    version_classes(X, m) :: Vector{Int}

The interchangeability partition of the `m` versions whose conflicts matrix is
`X`, as a class id per version — ids numbered `1, 2, …` in the version order of
each class's first member.

Per the manual's Theory section (lemma on interchangeable versions), two
versions of a package are *interchangeable* when they have the same dependency
set and every compatibility entry in the problem — their own and every other
package's — constrains them equally. In `PkgInfo` terms that is exactly row
equality, and the test is purely local to `X`: the dependency columns carry the
dependency set, and the interaction blocks carry every compatibility entry that
mentions the package, incoming as much as outgoing, because `pkg_info` records
each conflict symmetrically in both partners' matrices. (The incoming half is
not optional — another package's bound can admit one of two versions and not
the other, and then they are genuinely distinguishable.)

The partition is a property of the registry alone: it is independent of the
requirements, of the version ordering, and of any user constraints. That is
what puts it on the T1 side of the boundary, and what lets the collapse to
class representatives be a per-resolve step (see `class_representatives`).

The scan is word-parallel over columns: one pass computes, for every adjacent
pair of versions at once, whether their rows differ anywhere, which yields the
classes that are contiguous runs of versions — nearly all of them in practice,
since dependency sets and compat bounds vary in blocks. A second pass hashes
one row per run and merges runs whose rows compare equal, so the result is the
exact row-equality partition, hence invariant under reordering versions.
"""
function version_classes(X::BitMatrix, m::Int)
    version_classes!(Vector{Int}(undef, m), X, m)
end

function version_classes!(ids::Vector{Int}, X::BitMatrix, m::Int)
    resize!(ids, m)
    m > 0 || return ids
    n = size(X, 2) - 1
    W = col_words(X)
    chunks = X.chunks

    # pass 1: word-parallel adjacency differences. bit i-1 of `d` records
    # whether rows i and i+1 differ; within a word, row i's bit and row i+1's
    # bit are adjacent, so one shift-and-xor per word compares 64 pairs at
    # once (the topmost pair straddles into the next word)
    d = zeros(UInt64, W)
    @inbounds for k = 1:n
        base = (k - 1) * W
        for w = 1:W
            c = chunks[base + w]
            nxt = w < W ? chunks[base + w + 1] : UInt64(0)
            d[w] |= (c ⊻ (c >> 1)) ⊻ ((nxt & 1) << 63)
        end
    end
    runs = Int[1] # first version of each run of adjacent identical rows
    nc = 1
    ids[1] = 1
    @inbounds for i = 2:m
        if !iszero(d[((i - 2) >> 6) + 1] & (UInt64(1) << ((i - 2) & 63)))
            nc += 1
            push!(runs, i)
        end
        ids[i] = nc # = the index of i's run
    end
    nc > 1 || return ids

    # pass 2: merge runs whose rows are equal but not adjacent. hashes only
    # bucket the candidates; every merge is confirmed by comparing the rows
    buckets = Dict{UInt,Vector{Int}}()
    for (r, i) in enumerate(runs)
        push!(get!(() -> Int[], buckets, row_hash(X, i, n, W)), r)
    end
    canon = collect(1:nc) # each run's earliest run with an equal row
    merged = false
    for rs in values(buckets)
        length(rs) > 1 || continue
        sort!(rs)
        for a = 2:length(rs)
            for b = 1:a-1
                canon[rs[b]] == rs[b] || continue # not a class of its own
                rows_equal(X, runs[rs[b]], runs[rs[a]], n, W) || continue
                canon[rs[a]] = rs[b]
                merged = true
                break
            end
        end
    end
    merged || return ids

    # renumber the surviving classes in order of first appearance
    num = zeros(Int, nc)
    nc′ = 0
    for r = 1:nc
        canon[r] == r && (num[r] = (nc′ += 1))
    end
    @inbounds for i = 1:m
        ids[i] = num[canon[ids[i]]]
    end
    return ids
end

# renumber class ids to 1, 2, … in order of first appearance
function renumber_classes!(ids::Vector{Int})
    num = Dict{Int,Int}()
    for (i, c) in enumerate(ids)
        ids[i] = get!(() -> length(num) + 1, num, c)
    end
    return ids
end

function pkg_data(
    deps :: DepsProvider{P,D},
    reqs :: SetOrVec{P} = deps.packages;
) where {P,D}
    data = Dict{P,D}()
    work = Set(reqs)
    while !isempty(work)
        p = pop!(work)
        data_p = data[p] = pkg_data(deps, p)
        for (v, deps_pv) in data_p.depends, q in deps_pv
            q in keys(data) && continue
            push!(work, q)
        end
    end
    return data
end

function validate_pkg_data_consistency(data::AbstractDict{P,<:PkgData{P}}, reqs::SetOrVec{P}) where {P}
    available_packages = keys(data)

    # Check that all required packages exist
    for p in reqs
        if p ∉ available_packages
            throw(ArgumentError("Required package $p is not available in the package data"))
        end
    end

    # Check that all dependencies exist
    for (p, data_p) in data
        for deps_pv in values(data_p.depends)
            for q in deps_pv
                if q ∉ available_packages
                    throw(ArgumentError("Package $p depends on $q, but $q is not available in the package data"))
                end
            end
        end

        # Note: Compatibility constraints on unknown packages are allowed
    end
end

function pkg_info(
    deps :: DepsProvider{P},
    reqs :: SetOrVec{P} = deps.packages;
    filter :: Bool = true,
) where {P}
    data = pkg_data(deps, reqs)
    info = pkg_info(data, reqs; filter)
    return info
end

# a `Problem`'s requirements are all of it that reaches T1 — that is the point
# of the boundary — so these are pure conveniences for `pkg_info(_, prob.reqs)`.
# Elsewhere the `Problem` form is the primary one and the requirements form the
# convenience; here the dependence really does run the other way

pkg_info(
    deps :: DepsProvider{P},
    prob :: Problem{P};
    filter :: Bool = true,
) where {P} = pkg_info(deps, prob.reqs; filter)

pkg_info(
    data :: AbstractDict{P,<:PkgData{P}},
    prob :: Problem{P};
    filter :: Bool = true,
) where {P} = pkg_info(data, prob.reqs; filter)

"""
    pkg_info(data, reqs = all packages; filter = true) :: Dict{P,PkgInfo{P,V}}

The **T1 artifact**: the dependency closure of `reqs`, encoded as one `PkgInfo`
per package — conflict matrices, plus the
preprocessing that depends on the registry *alone*, namely the arc-consistency
prune (`mark_installable!`) and the interchangeability partition
(`version_classes`). Nothing here depends on the version ordering, on user
compat or pins, or (beyond the closure) on the requirements, so the result can
be computed once and shared by, or cached across, arbitrarily many resolves.

`filter = false` skips the arc-consistency prune, leaving the conflict matrices
exactly as the data spells them out.

Requirement-specific and constraint-specific shrinking — reachability and
redundancy elimination — happens per resolve instead, after the classes have
been refined and collapsed: see `prepare_pkg_info`.
"""
function pkg_info(
    data :: AbstractDict{P,<:PkgData{P,V}},
    reqs :: SetOrVec{P} = keys(data);
    filter :: Bool = true,
) where {P,V}
    # intern packages as integer indices once, up front: the build phases
    # below run entirely on vectors, and package names reappear only in
    # the final PkgInfo structures. ids follow name order (packages are
    # sorted), so id-sorted dependency and interaction lists coincide
    # with the name-sorted ones the output requires
    pkgs = sort!(collect(keys(data)))
    N = length(pkgs)
    ix = Dict{P,Int}(p => i for (i, p) in enumerate(pkgs))
    datas = [data[p] for p in pkgs]
    nver = Int[length(d.versions) for d in datas]

    # required packages must exist
    for r in reqs
        haskey(ix, r) || throw(ArgumentError(
            "Required package $r is not available in the package data"))
    end

    # conflict masks per distinct compat entry: for each partner q with a
    # conflicting version, the bitmask of conflicting versions of q
    MasksT = Vector{Tuple{Int,Vector{UInt64}}}
    function compute_masks(comp_pv)
        masks = MasksT()
        for (q, comp_pvq) in comp_pv
            qi = get(ix, q, 0)
            # compatibility constraints on unknown packages are allowed
            qi == 0 && continue
            vers = datas[qi].versions
            mask = zeros(UInt64, (length(vers) + 63) >> 6)
            found = false
            for (j, w) in enumerate(vers)
                w in comp_pvq && continue
                mask[((j - 1) >> 6) + 1] |= UInt64(1) << ((j - 1) & 63)
                found = true
            end
            found && push!(masks, (qi, mask))
        end
        return masks
    end

    # compute interactions between packages, keeping each version's
    # masks around for the bit-setting pass below. data providers share
    # compat dicts across a package's versions, so each package has only
    # a handful of distinct entries — a linear identity scan finds them
    # far more cheaply than hashing every version's dict
    vmasks = [Tuple{V,MasksT}[] for _ = 1:N]
    interacts = [Int[] for _ = 1:N]
    objs = Vector{Any}()
    mvals = Vector{MasksT}()
    for pi = 1:N
        interacts_p = interacts[pi]
        vm = vmasks[pi]
        empty!(objs)
        empty!(mvals)
        for (v, comp_pv) in datas[pi].compat
            idx = 0
            for ci = 1:length(objs)
                if objs[ci] === comp_pv
                    idx = ci
                    break
                end
            end
            if idx == 0
                push!(objs, comp_pv)
                push!(mvals, compute_masks(comp_pv))
                idx = length(objs)
            end
            masks = mvals[idx]
            isempty(masks) && continue
            push!(vm, (v, masks))
            for (qi, _) in masks
                (qi == pi || qi in interacts_p) && continue
                push!(interacts_p, qi)
                push!(interacts[qi], pi)
            end
        end
    end
    foreach(sort!, interacts)

    # dependency ids (validating that dependencies exist), interaction
    # block offsets, and the conflicts matrices
    depids = Vector{Vector{Int}}(undef, N)
    dobjs = Vector{Vector{Any}}(undef, N)       # distinct dep vectors
    dcols = Vector{Vector{Vector{Int}}}(undef, N) # their column indices
    offs = Vector{Vector{Int}}(undef, N)
    mats = Vector{BitMatrix}(undef, N)
    for pi = 1:N
        data_p = datas[pi]
        # union the versions' dependency sets; they are shared objects
        # (provider deduplication), so convert each distinct one once
        objs_p = Any[]
        for dv in values(data_p.depends)
            found = false
            for o in objs_p
                if o === dv
                    found = true
                    break
                end
            end
            found || push!(objs_p, dv)
        end
        ids = Int[]
        for dv in objs_p, q in dv
            qi = get(ix, q, 0)
            qi == 0 && throw(ArgumentError(
                "Package $(pkgs[pi]) depends on $q, but $q is not available in the package data"))
            push!(ids, qi)
        end
        ids = sort!(unique!(ids))
        depids[pi] = ids
        # per distinct dep vector, the dependency column indices
        # (self-dependencies get no bits, as before)
        cols_p = Vector{Int}[]
        for dv in objs_p
            cols = Int[]
            for q in dv
                qi = ix[q]
                qi == pi && continue
                push!(cols, searchsortedfirst(ids, qi))
            end
            push!(cols_p, cols)
        end
        dobjs[pi] = objs_p
        dcols[pi] = cols_p
        # interaction block offsets
        n = length(ids)
        off = Vector{Int}(undef, length(interacts[pi]))
        for (k, qi) in enumerate(interacts[pi])
            off[k] = n
            n += nver[qi]
        end
        offs[pi] = off
        # conflicts matrix: n+1 columns of padded_rows(m) bits each — rows
        # 1:m are versions, row m+1 holds column-active flags, the rest is
        # word-alignment padding (always zero); column n+1 holds the
        # version-active flags
        m = nver[pi]
        X = falses(padded_rows(m), n + 1)
        # mark all versions & columns as active; the column flags live at
        # a fixed bit of each column's chunk span, so set them directly
        # rather than through one strided BitArray write per column
        X[1:m, end] .= true
        W = col_words(X)
        wf = (m >> 6) + 1
        bf = UInt64(1) << (m & 63)
        ch = X.chunks
        @inbounds for j = 1:n
            ch[(j - 1) * W + wf] |= bf
        end
        mats[pi] = X
    end

    # block offset of package qi's versions in package pi's matrix
    block(pi::Int, qi::Int) =
        offs[pi][searchsortedfirst(interacts[pi], qi)]

    # initialize conflicts matrices
    for pi = 1:N
        data_p = datas[pi]
        X = mats[pi]
        V⁻¹ = Dict{V,Int}(v => i for (i, v) in enumerate(data_p.versions))
        # set dependency bits (column lists precomputed above)
        objs_p = dobjs[pi]
        cols_p = dcols[pi]
        for (v, deps_pv) in data_p.depends
            i = V⁻¹[v]
            idx = 0
            for ci = 1:length(objs_p)
                if objs_p[ci] === deps_pv
                    idx = ci
                    break
                end
            end
            for k in cols_p[idx]
                X[i, k] = true
            end
        end
        # set compatibility bits
        for (v, masks) in vmasks[pi]
            i = V⁻¹[v]
            Wx = col_words(X)
            wx = ((i - 1) >> 6) + 1
            bx = UInt64(1) << ((i - 1) & 63)
            xch = X.chunks
            for (qi, mask) in masks
                qi == pi && continue
                Y = mats[qi]
                b = block(pi, qi)
                c = block(qi, pi)
                # partner-side column is contiguous: blit the mask
                ybase = (c + i - 1) * col_words(Y)
                @inbounds for w = 1:length(mask)
                    Y.chunks[ybase + w] |= mask[w]
                end
                # own-side row: one bit per conflicting partner version,
                # written straight into the chunks (the row bit sits at a
                # fixed offset of each column's span)
                @inbounds for w = 1:length(mask)
                    z = mask[w]
                    while !iszero(z)
                        j = ((w - 1) << 6) + trailing_zeros(z) + 1
                        xch[(b + j - 1) * Wx + wx] |= bx
                        z &= z - 1
                    end
                end
            end
        end
    end

    # materialize the package-keyed PkgInfo structures. classes start out as
    # singletons (always a sound partition) and are computed for real below,
    # once the prune has settled which versions there are
    info = Dict{P,PkgInfo{P,V}}()
    for pi = 1:N
        D = pkgs[depids[pi]]
        T = Dict{P,Int}(
            pkgs[qi] => offs[pi][k]
            for (k, qi) in enumerate(interacts[pi]))
        info[pkgs[pi]] = PkgInfo{P,V}(
            datas[pi].versions, D, T, mats[pi], collect(1:nver[pi]))
    end

    # the T1 preprocessing: arc consistency, then interchangeability. Both are
    # properties of the registry alone — the arc-consistency test deletes
    # versions one of whose dependencies has no versions at all, which is a
    # sound approximation of deleting model-free versions and, per the manual's
    # Theory section, exactly the kind of deletion that is safe for *every*
    # ordering and requirement set; the classes are row equality. Reachability
    # and redundancy elimination are not of that kind — both are stated in
    # terms of the version ordering, and reachability in terms of the
    # requirements too — so they now run per resolve, in `prepare_pkg_info`.
    if filter
        mark_installable!(info)
        drop_unmarked!(info)
    end
    for info_p in values(info)
        version_classes!(info_p.classes, info_p.conflicts,
                         length(info_p.versions))
    end

    return info
end
