# Type Parameters:
#  P = package type (String)
#  V = version type (VersionNumber)
#  S = version set type (VersionSpec)

struct PkgInfo{P,V,S}
    versions :: Vector{V}
    depends  :: Dict{V, Vector{P}}
    compat   :: Dict{V, Dict{P, S}}
end

struct DepsProvider{P,V,S,F<:Function}
    provider :: F
end

DepsProvider{P,V,S}(provider::Function) where {P,V,S} =
    DepsProvider{P,V,S,typeof(provider)}(provider)

(deps::DepsProvider{P,V,S,F})(pkg::P) where {P,V,S,F<:Function} =
    deps.provider(pkg) :: PkgInfo{P,V,S}

# find all packages and versions that might be needed for resolution

function find_packages(
    deps :: DepsProvider{P,V,S},
    reqs :: Vector{P},
) where {P,V,S}
    pkgs = Dict{P, PkgInfo{P,V,S}}()
    work = Set(reqs)
    while !isempty(work)
        pkg = pop!(work)
        @assert pkg ∉ keys(pkgs)
        info = pkgs[pkg] = deps(pkg)
        for pkgs′ in values(info.depends), pkg′ in pkgs′
            pkg′ in keys(pkgs) && continue
            push!(work, pkg′)
        end
    end
    return pkgs
end

# filter down to versions that resolve might actually pick

"""
    find_reachable(deps, reqs) :: Dict{P, Int}

This function finds a minimal "reachable" subset of package and versions that
could appear in pareto-optimal solutions to version resolution for the given set
of required "root" packages, using the following recursive logic:

- P in reqs => P[1] reachable
- P[i] reachable & P[i] depends on D => D[1] reachable
- P[i] reachable & P[i] conflicts w/ something reachable => P[i+1] reachable

The function returns a dictionary mapping packages to the maximum version index
of that package that could be reached in an optimal solution. If a pacakge
cannot appear in an optimal solution, it will not appear in this dictionary.
"""
function find_reachable(
    pkgs :: Dict{P, PkgInfo{P,V,S}},
    reqs :: Vector{P},
) where {P,V,S}
    reach = Dict{P,Int}(p => 0 for p in reqs)
    queue = Dict{P,Int}(p => 1 for p in reqs)

    function add_queue!(dep::P, k::Int)
        i = get(reach, dep, 0)
        i ≥ k && return false
        j = get(queue, dep, 0)
        j ≥ k && return false
        queue[dep] = k
        return true
    end

    while !isempty(queue)
        # get unprocessed package + version
        (pkg, k) = pop!(queue)
        pkg_info = pkgs[pkg]
        j = get(reach, pkg, 0)
        k = min(k, length(pkg_info.versions))
        reach[pkg] = k
        for i = j+1:k
            pkg_v = pkg_info.versions[i]
            # look at dependencies
            pkg_v in keys(pkg_info.depends) &&
            for dep in pkg_info.depends[pkg_v]
                add_queue!(dep, 1)
            end
            # look at conflicts
            pkg_v in keys(pkg_info.compat) &&
            for (dep, spec) in pkg_info.compat[pkg_v]
                dep_info = pkgs[dep]
                l = get(reach, dep, 0)
                for (m, dep_v) in enumerate(dep_info.versions)
                    m ≤ l || break
                    # check compatibilty
                    compatible = dep_v in spec
                    if dep_v in keys(dep_info.compat)
                        dep_v_compat = dep_info.compat[dep_v]
                        if pkg in keys(dep_v_compat)
                            compatible &= pkg_v in dep_v_compat[pkg]
                        end
                    end
                    compatible && continue
                    # incompatible
                    add_queue!(pkg, i+1)
                    add_queue!(dep, m+1)
                end
            end
        end
    end

    return reach
end

function filter_reachable!(
    pkgs  :: Dict{P, PkgInfo{P,V,S}},
    reach :: Dict{P, Int},
) where {P,V,S}
    for pkg in sort!(collect(keys(reach)))
        k = reach[pkg]
        info = pkgs[pkg]
        vers = info.versions
        @assert k ≤ length(vers)
        for j = k+1:length(vers)
            ver = vers[j]
            delete!(info.depends, ver)
            delete!(info.compat, ver)
        end
        resize!(vers, k)
    end
    filter!(pkgs) do (pkg, info)
        pkg in keys(reach) && !isempty(info.versions)
    end
end

function filter_reachable!(
    pkgs :: Dict{P, PkgInfo{P,V,S}},
    reqs :: Vector{P},
) where {P,V,S}
    reach = find_reachable(pkgs, reqs)
    filter_reachable!(pkgs, reach)
end

# two packages interact if there is some conflict between them

function find_interacts(
    pkgs :: Dict{P, PkgInfo{P,V,S}},
) where {P,V,S}
    interacts = Dict{P,Vector{P}}(p => P[] for p in keys(pkgs))
    for (pkg₁, info₁) in pkgs
        interact₁ = interacts[pkg₁]
        for ver₁ in info₁.versions
            ver₁ in keys(info₁.compat) || continue
            compat₁ = info₁.compat[ver₁]
            for (pkg₂, spec₁) in compat₁
                pkg₂ in interact₁ && continue
                interact₂ = interacts[pkg₂]
                for ver₂ in pkgs[pkg₂].versions
                    if ver₂ ∉ spec₁
                        push!(interact₁, pkg₂)
                        push!(interact₂, pkg₁)
                        break
                    else
                        compat₂ = pkgs[pkg₂].compat
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
    filter!(interacts) do (pkg, ix)
        !isempty(ix)
    end
    foreach(sort!, values(interacts))
    return interacts
end

# compute the set of conflicts between package versions

struct Conflicts{P}
    depends   :: Vector{P}
    interacts :: Dict{P, Int}
    conflicts :: Matrix{Bool}
end

function find_conflicts(
    pkgs :: Dict{P, PkgInfo{P,V,S}},
    active :: Bool,
    p :: P,
    intx :: Vector{P},
) where {P,V,S}
    # look up some stuff about p
    info = pkgs[p]
    vers_p = info.versions
    comp_p = info.compat
    # collect all dependency packages
    dx = P[]
    for (v, deps) in info.depends
        union!(dx, deps)
    end
    sort!(dx)
    ix = Dict{P, Int}(p => 0 for p in intx)
    # compute interactions matrix (m × n)
    m = length(vers_p) # per package version
    n = length(dx) + # per dependency package + interacting package version
        sum(init=0, length(pkgs[q].versions) for q in intx)
    X = fill(false, m + active, n + active) # conflicts & actives
    for (i, v) in enumerate(vers_p)
        deps_pv = info.depends[v]
        comp_pv = info.compat[v]
        for (j, q) in enumerate(dx)
            X[i, j] = q ∈ deps_pv
        end
        b = length(dx)
        for q in intx
            ix[q] = b
            vers_q = pkgs[q].versions
            comp_q = pkgs[q].compat
            for (j, u) in enumerate(vers_q)
                comp_qu = comp_q[u]
                X[i, b + j] =
                    v ∈ keys(comp_p) &&
                    q ∈ keys(comp_pv) &&
                    u ∉ comp_pv[q] ||
                    u ∈ keys(comp_q) &&
                    p ∈ keys(comp_qu) &&
                    v ∉ comp_qu[p]
            end
            b += length(vers_q)
        end
    end
    # initially all versions are active
    if active
        X[1:m, end] .= true
        # columns can start inactive if they have no conflicts
        for j = 1:n
            X[end, j] = any(X[i, j] for i = 1:m)
        end
    end
    return Conflicts(dx, ix, X)
end

function find_conflicts(
    pkgs :: Dict{P,PkgInfo{P,V,S}},
    active :: Bool = false,
) where {P,V,S}
    interacts = find_interacts(pkgs)
    ∅ = P[] # reuse empty vector
    Dict{P, Conflicts{P}}(
        p => find_conflicts(pkgs, active, p, get(interacts, p, ∅))
        for p in keys(pkgs)
    )
end

function filter_redundant!(
    pkgs :: Dict{P,PkgInfo{P,V,S}},
    cx   :: Dict{P,Conflicts{P}} = find_conflicts(pkgs, true),
) where {P,V,S}
    work = copy(keys(cx))
    names = sort!(collect(work))
    # some work vectors
    J = Int[] # active versions vector
    K = Int[] # active conflicts vector
    R = Int[] # redundant indices vector
    sort!(names, by = p -> length(cx[p].interacts))
    for p in Iterators.cycle(names)
        isempty(work) && break
        p in work || continue
        delete!(work, p)
        # get conflicts & dimensions
        X = cx[p].conflicts
        m = size(X, 1) - 1
        m > 1 || continue # unique version cannot be reundant
        n = size(X, 2) - 1
        # active indices
        append!(empty!(J), (j for j = 1:m if X[j, end]))
        length(J) > 1 || continue # unique version cannot be reundant
        append!(empty!(K), (k for k = 1:n if X[end, k]))
        # find redundant versions
        empty!(R)
        for j in J
            for i in J
                i < j || break
                i ∈ R && continue
                all(!X[i, k] | X[j, k] for k in K) || continue
                # an earlier version is strictly more compatible
                # i.e. i < j and X[i, k] => X[j, k] for all k
                # therefore i will always be chosen instead of j
                push!(R, j)
                break
            end
        end
        isempty(R) && continue
        # deactivate redundant versions
        X[R, end] .= false
        for q in keys(cx[p].interacts)
            b = cx[q].interacts[p]
            cx[q].conflicts[end, b .+ R] .= false
            push!(work, q) # can create new redundancies
        end
    end
    for (p, info) in pkgs
        X = cx[p].conflicts
        m = length(info.versions)
        append!(empty!(R), (j for j = 1:m if !X[j, end]))
        for (i, v) in enumerate(info.versions)
            i in R || continue
            delete!(info.depends, v)
            delete!(info.compat, v)
        end
        deleteat!(info.versions, R)
    end
end
