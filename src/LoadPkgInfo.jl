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
    data = Dict{P, PkgData{P,V,S}}()
    work = Set(reqs)
    while !isempty(work)
        pkg = pop!(work)
        @assert pkg ∉ keys(data)
        entry = data[pkg] = deps(pkg)
        for entries′ in values(entry.depends), pkg′ in entries′
            pkg′ in keys(data) && continue
            push!(work, pkg′)
        end
    end
    # now, compute interactions between packages
    interacts = Dict{P,Vector{P}}(p => P[] for p in keys(data))
    for (pkg₁, data₁) in data
        interact₁ = interacts[pkg₁]
        for ver₁ in data₁.versions
            ver₁ in keys(data₁.compat) || continue
            compat₁ = data₁.compat[ver₁]
            for (pkg₂, spec₁) in compat₁
                pkg₂ in interact₁ && continue
                interact₂ = interacts[pkg₂]
                for ver₂ in data[pkg₂].versions
                    if ver₂ ∉ spec₁
                        push!(interact₁, pkg₂)
                        push!(interact₂, pkg₁)
                        break
                    else
                        compat₂ = data[pkg₂].compat
                        pkg₁ in keys(compat₂) || continue
                        spec₂ = compat₂[pkg₁]
                        if ver₁ ∉ spec₂
                            push!(interact₁, pkg₂)
                            push!(interact₂, pkg₁)
                            break
                        end
                    end
                end
            end
        end
    end
    foreach(sort!, values(interacts))
    # construct dict of PkgInfo structs
    info = Dict{P, PkgInfo{P,V}}()
    for (p, entry) in data
        interacts_p = interacts[p]
        vers_p = entry.versions
        deps_p = entry.depends
        comp_p = entry.compat
        # collect all dependency packages
        dx = P[]
        for (v, deps) in deps_p
            union!(dx, deps)
        end
        sort!(dx)
        ix = Dict{P, Int}(p => 0 for p in interacts_p)
        # compute interactions matrix (m × n)
        m = length(vers_p) # per package version
        n = length(dx) + # per dependency package + interacting package version
            sum(init=0, length(data[q].versions) for q in interacts_p)
        X = falses(m + 1, n + 1) # conflicts & active flags
        # empty deps & compat
        deps_∅ = Vector{P}()
        comp_∅ = Dict{P,S}()
        # main work loop
        for (i, v) in enumerate(vers_p)
            deps_pv = get(deps_p, v, deps_∅)
            comp_pv = get(comp_p, v, comp_∅)
            # dependencies
            for (j, q) in enumerate(dx)
                X[i, j] = q ∈ deps_pv
            end
            # conflicts
            b = length(dx)
            for q in interacts_p
                ix[q] = b
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
        info[p] = PkgInfo(vers_p, dx, ix, X)
    end
    if filter
        mark_reachable!(info, reqs)
        mark_necessary!(info)
        drop_unmarked!(info)
    end
    return info
end
