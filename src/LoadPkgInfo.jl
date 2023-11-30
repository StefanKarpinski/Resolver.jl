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

function load_pkg_info(
    deps :: DepsProvider{P,V,S},
    reqs :: SetOrVec{P} = deps.packages;
    filter :: Bool = true,
) where {P,V,S}
    # first, load dict of PkgData structs
    data = Dict{P,PkgData{P,V,S}}()
    work = Set(reqs)
    while !isempty(work)
        p = pop!(work)
        @assert p ∉ keys(data)
        data_p = data[p] = deps(p)
        for (v, deps_pv) in data_p.depends, q in deps_pv
            q in keys(data) && continue
            push!(work, q)
        end
    end
    # now, compute interactions between packages
    interacts = Dict{P,Vector{P}}(p => P[] for p in keys(data))
    for (p, data_p) in data
        interacts_p = interacts[p]
        for v in data_p.versions
            v in keys(data_p.compat) || continue
            compat_pv = data_p.compat[v]
            for (q, spec_q) in compat_pv
                q in interacts_p && continue
                interacts_q = interacts[q]
                data_q = data[q]
                for w in data_q.versions
                    if w ∉ spec_q
                        push!(interacts_p, q)
                        push!(interacts_q, p)
                        break
                    else
                        compat_q = data_q.compat
                        p in keys(compat_q) || continue
                        if v ∉ compat_q[p]
                            push!(interacts_p, q)
                            push!(interacts_q, p)
                            break
                        end
                    end
                end
            end
        end
    end
    foreach(sort!, values(interacts))
    # construct dict of PkgInfo structs
    info = Dict{P,PkgInfo{P,V}}()
    for (p, data_p) in data
        vers_p = data_p.versions
        deps_p = data_p.depends
        comp_p = data_p.compat
        # collect dependencies across all versions
        deps_pa = sort!(reduce(union!, values(deps_p), init=P[]))
        interacts_p = Dict{P, Int}(p => 0 for p in interacts[p])
        # compute interactions matrix (m × n)
        m = length(vers_p) # per package version
        n = length(deps_pa) + # per dependency + per interacting version
            sum(init=0, length(data[q].versions) for q in keys(interacts_p))
        X = falses(m + 1, n + 1) # conflicts & active flags
        # empty deps & compat
        deps_∅ = Vector{P}()
        comp_∅ = Dict{P,S}()
        # main work loop
        for (i, v) in enumerate(vers_p)
            deps_pv = get(deps_p, v, deps_∅)
            comp_pv = get(comp_p, v, comp_∅)
            # dependencies
            for (j, q) in enumerate(deps_pa)
                X[i, j] = q ∈ deps_pv
            end
            # conflicts
            b = length(deps_pa)
            for q in interacts[p]
                interacts_p[q] = b
                info_q = data[q]
                vers_q = info_q.versions
                comp_q = info_q.compat
                for (j, w) in enumerate(vers_q)
                    comp_qw = get(comp_q, w, comp_∅)
                    X[i, b + j] =
                        v ∈ keys(comp_p) &&
                        q ∈ keys(comp_pv) &&
                        w ∉ comp_pv[q] ||
                        w ∈ keys(comp_q) &&
                        p ∈ keys(comp_qw) &&
                        v ∉ comp_qw[p]
                end
                b += length(vers_q)
            end
        end
        # mark all versions as active
        X[1:m, end] .= true
        X[end, 1:n] .= true
        X[end, end] = false
        # add the PkgInfo struct to dict
        info[p] = PkgInfo(vers_p, deps_pa, interacts_p, X)
    end
    filter && filter_pkg_info!(info, reqs)
    return info
end
