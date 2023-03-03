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
    pkgs = Dict{P,PkgInfo{P,V,S}}()
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
    pkgs :: Dict{P,PkgInfo{P,V,S}},
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
    pkgs::Dict{P,PkgInfo{P,V,S}},
    reach::Dict{P,Int},
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
end

# eliminate redundant versions

"""
    filter_redundant!(pkgs)

"""

# compute the set of conflicts between package versions

struct Conflicts{P,V}
    versions  :: Vector{Pair{P,V}}
    conflicts :: Vector{Vector{Int}}
end

function find_conflicts(pkgs::Dict{P,PkgInfo{P,V,S}}) where {P,V,S}
    conflicts = Dict{Tuple{P,P}, Vector{Tuple{V,V}}}()

    function add_conflict!(p₁::P, v₁::V, p₂::P, v₂::V)
        p₂ < p₁ && return add_conflict!(p₂, v₂, p₁, v₁)
        pairs = get!(()->valtype(conflicts)(), conflicts, (p₁, p₂))
        (v₁, v₂) in pairs || push!(pairs, (v₁, v₂))
    end

    for (pkg₁, info₁) in sort!(collect(pkgs), by=first),
        ver₁ in info₁.versions
        ver₁ in keys(info₁.compat) || continue
        compat₁ = info₁.compat[ver₁]
        for (pkg₂, spec₁) in compat₁
            @show pkg₁, pkg₂
            for ver₂ in pkgs[pkg₂].versions
                if ver₂ ∉ spec₁
                    add_conflict!(pkg₁, ver₁, pkg₂, ver₂)
                else
                    compat₂ = pkgs[pkg₂].compat
                    pkg₁ in keys(compat₂) || continue
                    spec₂ = compat₂[pkg₁]
                    if ver₁ ∉ spec₂
                        add_conflict!(pkg₁, ver₁, pkg₂, ver₂)
                    end
                end
            end
        end
        @show length(conflicts)
        @show sum(length, values(conflicts))
    end
    filter!(conflicts) do (pkgs, vers)
        !isempty(vers)
    end
    return conflicts
end
