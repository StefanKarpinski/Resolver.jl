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
    deps :: DepsProvider{P,D},
    reqs :: SetOrVec{P} = deps.packages;
    filter :: Bool = true,
) where {P,D}
    # first, load dict of PkgData structs
    data = Dict{P,D}()
    work = Set(reqs)
    @timeit "load pkg data" while !isempty(work)
        p = pop!(work)
        @assert p ∉ keys(data)
        data_p = data[p] = pkg_data(deps, p)
        for (v, deps_pv) in data_p.depends, q in deps_pv
            q in keys(data) && continue
            push!(work, q)
        end
    end
    return make_pkg_info(data, reqs; filter)
end

function make_pkg_info(
    data :: AbstractDict{P, <:PkgData{P,V}},
    reqs :: SetOrVec{P} = keys(data);
    filter :: Bool = true,
) where {P,V}
    # compute interactions between packages
    interacts = Dict{P,Vector{P}}(p => P[] for p in keys(data))
    @timeit "compute interacts" for (p, data_p) in data
        interacts_p = interacts[p]
        for (v, comp_pv) in data_p.compat
            # don't combine loops--it changes what continue does
            for (q, comp_pvq) in comp_pv
                (q == p || q in interacts_p) && continue
                interacts_q = interacts[q]
                for w in data[q].versions
                    w in comp_pvq && continue
                    push!(interacts_p, q)
                    push!(interacts_q, p)
                    break
                end
            end
        end
    end
    foreach(sort!, values(interacts))

    # construct dict of PkgInfo structs
    info = Dict{P,PkgInfo{P,V}}()
    @timeit "construct pkg info" for (p, data_p) in data
        D = sort!(reduce(union!, values(data_p.depends), init=P[]))
        T = Dict{P,Int}(p => 0 for p in interacts[p])
        n = length(D)
        for q in interacts[p]
            T[q] = n
            n += length(data[q].versions)
        end
        # conflicts matrix (m + 1) × (n + 1)
        m = length(data_p.versions)
        X = falses(m + 1, n + 1)
        # mark all versions as active
        X[1:m, end] .= true
        X[end, 1:n] .= true
        # add the PkgInfo struct to dict
        info[p] = PkgInfo(data_p.versions, D, T, X)
    end

    # initialize conflicts matrices
    @timeit "initialize conflicts" for (p, info_p) in info
        X = info_p.conflicts
        V⁻¹ = Dict{V,Int}(v => i for (i, v) in enumerate(info_p.versions))
        D⁻¹ = Dict{P,Int}(q => j for (j, q) in enumerate(info_p.depends))
        data_p = data[p]
        # set dependency bits
        for (v, deps_pv) in data_p.depends
            i = V⁻¹[v]
            for q in deps_pv
                X[i, D⁻¹[q]] = true
            end
        end
        # set compatibility bits
        for (v, comp_pv) in data_p.compat
            i = V⁻¹[v]
            for (q, comp_pvq) in comp_pv
                haskey(info_p.interacts, q) || continue
                info_q = info[q]
                Y = info_q.conflicts
                b = info_p.interacts[q]
                c = info_q.interacts[p]
                for (j, w) in enumerate(info_q.versions)
                    w in comp_pvq && continue
                    X[i, b + j] = true
                    Y[j, c + i] = true
                end
            end
        end
    end

    # only keep reachable, necessary packages & versions
    filter && @timeit "filter pkg info" filter_pkg_info!(info, reqs)

    return info
end
