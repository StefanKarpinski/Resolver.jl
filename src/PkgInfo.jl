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
    validate_pkg_data_consistency(data, reqs)
    # conflict masks per distinct compat entry: for each partner q with a
    # conflicting version, the bitmask of conflicting versions of q. data
    # providers commonly share compat dicts across versions, so each
    # distinct entry is scanned against the partners' versions only once
    memo = IdDict{Any, Dict{P, Vector{UInt64}}}()
    function conflict_masks(comp_pv)
        get!(memo, comp_pv) do
            masks = Dict{P, Vector{UInt64}}()
            for (q, comp_pvq) in comp_pv
                haskey(data, q) || continue # unreachable (weak dep)
                vers = data[q].versions
                mask = zeros(UInt64, (length(vers) + 63) >> 6)
                found = false
                for (j, w) in enumerate(vers)
                    w in comp_pvq && continue
                    mask[((j - 1) >> 6) + 1] |= UInt64(1) << ((j - 1) & 63)
                    found = true
                end
                found && (masks[q] = mask)
            end
            masks
        end
    end

    # compute interactions between packages
    interacts = Dict{P,Vector{P}}(p => P[] for p in keys(data))
    for (p, data_p) in data
        interacts_p = interacts[p]
        for (v, comp_pv) in data_p.compat
            for q in keys(conflict_masks(comp_pv))
                (q == p || q in interacts_p) && continue
                push!(interacts_p, q)
                push!(interacts[q], p)
            end
        end
    end
    foreach(sort!, values(interacts))

    # construct dict of PkgInfo structs
    info = Dict{P,PkgInfo{P,V}}()
    for (p, data_p) in data
        D = sort!(reduce(union!, values(data_p.depends), init=P[]))
        T = Dict{P,Int}(p => 0 for p in interacts[p])
        n = length(D)
        for q in interacts[p]
            T[q] = n
            n += length(data[q].versions)
        end
        # conflicts matrix: n+1 columns of padded_rows(m) bits each — rows
        # 1:m are versions, row m+1 holds column-active flags, the rest is
        # word-alignment padding (always zero); column n+1 holds the
        # version-active flags
        m = length(data_p.versions)
        X = falses(padded_rows(m), n + 1)
        # mark all versions & columns as active
        X[1:m, end] .= true
        X[m+1, 1:n] .= true
        # add the PkgInfo struct to dict
        info[p] = PkgInfo{P,V}(data_p.versions, D, T, X)
    end

    # initialize conflicts matrices
    for (p, info_p) in info
        X = info_p.conflicts
        V⁻¹ = Dict{V,Int}(v => i for (i, v) in enumerate(info_p.versions))
        D⁻¹ = Dict{P,Int}(q => j for (j, q) in enumerate(info_p.depends))
        data_p = data[p]
        # set dependency bits
        for (v, deps_pv) in data_p.depends
            i = V⁻¹[v]
            for q in deps_pv
                q == p && continue
                X[i, D⁻¹[q]] = true
            end
        end
        # set compatibility bits
        for (v, comp_pv) in data_p.compat
            i = V⁻¹[v]
            for (q, mask) in conflict_masks(comp_pv)
                q == p && continue
                info_q = info[q]
                Y = info_q.conflicts
                b = info_p.interacts[q]
                c = info_q.interacts[p]
                # partner-side column is contiguous: blit the mask
                ybase = (c + i - 1) * col_words(Y)
                @inbounds for w = 1:length(mask)
                    Y.chunks[ybase + w] |= mask[w]
                end
                # own-side row: one bit per conflicting partner version
                @inbounds for w = 1:length(mask)
                    z = mask[w]
                    while !iszero(z)
                        X[i, b + ((w - 1) << 6) + trailing_zeros(z) + 1] = true
                        z &= z - 1
                    end
                end
            end
        end
    end

    # only keep reachable, necessary packages & versions
    filter && filter_pkg_info!(info, reqs)

    return info
end
