struct PkgInfo{P,V}
    versions  :: Vector{V}
    depends   :: Vector{P}
    interacts :: Dict{P, Int}
    conflicts :: BitMatrix
end

function Base.:(==)(a::PkgInfo, b::PkgInfo)
    a.versions  == b.versions  &&
    a.depends   == b.depends   &&
    a.interacts == b.interacts &&
    a.conflicts == b.conflicts
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

    # materialize the package-keyed PkgInfo structures
    info = Dict{P,PkgInfo{P,V}}()
    for pi = 1:N
        D = pkgs[depids[pi]]
        T = Dict{P,Int}(
            pkgs[qi] => offs[pi][k]
            for (k, qi) in enumerate(interacts[pi]))
        info[pkgs[pi]] = PkgInfo{P,V}(datas[pi].versions, D, T, mats[pi])
    end

    # only keep reachable, necessary packages & versions
    filter && filter_pkg_info!(info, reqs)

    return info
end
