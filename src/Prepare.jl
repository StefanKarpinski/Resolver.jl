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
    pkgs  :: Dict{P,PkgInfo{P,V,S}},
    reach :: Dict{P,Int},
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

# two packages interact if there is some conflict between them

function find_interactions(
    pkgs :: Dict{P,PkgInfo{P,V,S}},
) where {P,V,S}
    interacts = Dict{P,Vector{P}}(
        p => Vector{P}() for p in keys(pkgs)
    )
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

function find_conflicts(
    pkgs      :: Dict{P,PkgInfo{P,V,S}},
    interacts :: Dict{P,Vector{P}},
) where {P,V,S}
    conflicts = Dict{P, BitMatrix}()

    for p₁ in sort!(collect(keys(interacts)))
        vers₁ = pkgs[p₁].versions
        comp₁ = pkgs[p₁].compat
        m = length(vers₁)
        N = sum(length(pkgs[p₂].versions) for p₂ in interacts[p₁])
        conflicts[p₁] = valtype(conflicts)(undef, m, N)
        b = 0
        for p₂ in interacts[p₁]
            vers₂ = pkgs[p₂].versions
            comp₂ = pkgs[p₂].compat
            n = length(vers₂)
            conflicts[p₁][:, b .+ (1:n)] = [
                v₁ ∈ keys(comp₁) &&
                p₂ ∈ keys(comp₁[v₁]) &&
                v₂ ∉ comp₁[v₁][p₂] ||
                v₂ ∈ keys(comp₂) &&
                p₁ ∈ keys(comp₂[v₂]) &&
                v₁ ∉ comp₂[v₂][p₁] for
                (i₁, v₁) in enumerate(vers₁),
                (i₂, v₂) in enumerate(vers₂)
            ]
            b += n
        end
    end

    return conflicts
end

function find_redundant(
    conflicts :: Dict{P,<:AbstractMatrix{Bool}},
    dirty     :: AbstractSet{P} = keys(conflicts),
) where P
    redundant = Dict{P,Vector{Int}}()
    for pkg in dirty
        X = conflicts[pkg]
        R = Int[]
        for j = 2:size(X, 1)
            for i = 1:j-1
                i in R && continue
                if all(!X[i, k] | X[j, k] for k = 1:size(X, 2))
                    # an earlier version is strictly more compatible
                    # i.e. i < j and X[i, k] => X[j, k] for all k
                    push!(R, j)
                    break
                end
            end
        end
        isempty(R) && continue
        redundant[pkg] = R
    end
    return redundant
end

struct ConflictInfo{P}
    interacts :: Vector{P}
    conflicts :: BitMatrix
end

function filter_redundant!(
    pkgs :: Dict{P,PkgInfo{P,V,S}},
) where {P,V,S}
    while true
        interacts = find_interactions(pkgs)
        conflicts = find_conflicts(pkgs, interacts)
        redundant = find_redundant(conflicts)
        isempty(redundant) && return Dict{P,ConflictInfo{P}}(
            p => ConflictInfo(ix, conflicts[p]) for (p, ix) in interacts
        )
        while isempty(redundant)
            dirty = filter_redundant!(pkgs, conflicts, redundant)
            redundant = find_redundant(conflicts, dirty)
        end
    end
end

# function filter_redundant!(
#     pkgs      :: Dict{P,PkgInfo{P,V,S}},
#     conflicts :: Dict{P,<:AbstractMatrix{Bool}},
# ) where {P,V,S}
#     dirty = copy(keys(conflicts))
#     while !isempty(dirty)
#         redundant = find_redundant(conflicts, dirty)
#         dirty = filter_redundant!(pkgs, conflicts, redundant)
#     end
# end

function filter_redundant!(
    pkgs      :: Dict{P,PkgInfo{P,V,S}},
    conflicts :: Dict{P,<:AbstractMatrix{Bool}},
    redundant :: Dict{P,Vector{Int}},
) where {P,V,S}
    dirty = Set{P}()
    # first filter conflicts
    for (p₁, R₁) in redundant
        N = 0
        R₂ = Int[]
        for p₂ in interacts[p₁]
            push!(dirty, p₂)
            if p₂ in keys(redundant)
                append!(R₂, N .+ redundant[p₂])
            end
            N += length(pkgs[p₂].versions)
        end
        m = length(pkgs[p₁].versions)
        I = setdiff(1:m, R₁)
        J = setdiff(1:N, R₂)
        conflicts[p₁] = conflicts[p₁][I,J]
    end
    # next filter pkgs to match
    for (p, R) in redundant
        info = pkgs[p]
        for (i, v) in enumerate(info.versions)
            i in R || continue
            delete!(info.depends, v)
            delete!(info.compat, v)
        end
        deleteat!(info.versions, R)
        @assert !isempty(info.versions)
    end
    # return potentially affected packages
    return dirty
end
