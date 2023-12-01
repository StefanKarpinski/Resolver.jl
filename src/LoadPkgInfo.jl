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

using TimerOutputs

function load_pkg_info(
    deps :: DepsProvider{P,V,S},
    reqs :: SetOrVec{P} = deps.packages;
    filter :: Bool = true,
) where {P,V,S}
    # first, load dict of PkgData structs
    data = Dict{P,PkgData{P,V,S}}()
    work = Set(reqs)
    @timeit "load pkg data" while !isempty(work)
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
    # POSSIBLE TODO: more efficient data structure
    #   - sort names, use Vector{Vector{Int}}
    #   - indices imply names
    #   - 50% less memory
    #   - sorted construction
    @timeit "compute interacts" for (p, data_p) in data
        interacts_p = interacts[p]
        for v in data_p.versions
            v in keys(data_p.compat) || continue
            for (q, comp_pvq) in data_p.compat[v]
                q in interacts_p && continue
                interacts_q = interacts[q]
                data_q = data[q]
                for w in data_q.versions
                    if w ∉ comp_pvq
                        push!(interacts_p, q)
                        push!(interacts_q, p)
                        break
                    else
                        comp_qw = data_q.compat[w]
                        p in keys(comp_qw) || continue
                        if v ∉ comp_qw[p]
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
    @timeit "construct PkgInfo dict" for (p, data_p) in data
        vers_p = data_p.versions
        deps_p = data_p.depends
        comp_p = data_p.compat
        # collect dependencies across all versions
        deps_pa = sort!(reduce(union!, values(deps_p), init=P[]))
        deps_pd = Dict{P,Int}(p => i for (i, p) in enumerate(deps_pa))
        interacts_p = Dict{P, Int}(p => 0 for p in interacts[p])
        # compute interactions matrix (m × n)
        m = length(vers_p) # per package version
        n = length(deps_pd) + # per dependency + per interacting version
            sum(init=0, length(data[q].versions) for q in keys(interacts_p))
        X = falses(m + 1, n + 1) # conflicts & active flags
        # set dependency bits
        @timeit "deps" for (i, v) in enumerate(vers_p)
            for q in deps_p[v]
                j = deps_pd[q]
                X[i, j] = true
            end
        end
        # empty deps & compat
        comp_∅ = Dict{P,S}()
        # set conflict bits
        b = length(deps_pd)
        @timeit "conflicts" for q in interacts[p]
            info_q = data[q]
            vers_q = info_q.versions
            comp_q = info_q.compat
            for (j, w) in enumerate(vers_q)
                comp_qw = get(comp_q, w, comp_∅)
                for (i, v) in enumerate(vers_p)
                    comp_pv = get(comp_p, v, comp_∅)
                    X[i, b + j] =
                        v ∈ keys(comp_p) &&
                        q ∈ keys(comp_pv) &&
                        w ∉ comp_pv[q] ||
                        w ∈ keys(comp_q) &&
                        p ∈ keys(comp_qw) &&
                        v ∉ comp_qw[p]
                end
            end
            interacts_p[q] = b
            b += length(vers_q)
        end
        # mark all versions as active
        X[1:m, end] .= true
        X[end, 1:n] .= true
        X[end, end] = false
        # add the PkgInfo struct to dict
        info[p] = PkgInfo(vers_p, deps_pa, interacts_p, X)
    end
    @timeit "filter pkg info" filter && filter_pkg_info!(info, reqs)
    return info
end
