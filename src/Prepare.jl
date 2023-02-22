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

# compute the set of conflicts between package versions

function find_conflicts(
    pkgs :: Dict{P,PkgInfo{P,V,S}},
) where {P,V,S}
    conflicts = Dict{Tuple{P,P}, Set{Tuple{V,V}}}()
    add_conflict!(p₁, v₁, p₂, v₂) = p₂ < p₁ ? add_conflict!(p₂, v₂, p₁, v₁) :
        push!(get!(()->valtype(conflicts)(), conflicts, (p₁, p₂)), (v₁, v₂))

    for (pkg₁, info₁) in pkgs, ver₁ in info₁.versions
        ver₁ in keys(info₁.compat) || continue
        compat₁ = info₁.compat[ver₁]
        for (pkg₂, spec₁) in compat₁, ver₂ in pkgs[pkg₂].versions
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
    filter!(conflicts) do (pkgs, vers)
        !isempty(vers)
    end
    return conflicts
end

# filter down to versions that resolve might actually pick

struct PkgVer{P,V}
    version   :: V
    depends   :: Vector{P}
    conflicts :: Vector{Vector{V}}
end

"""
    reachable(deps::DepsProvider{P,V,S}, reqs::Vector{P})
        :: Dict{P, Vector{PkgVer{P,V}}

This function finds a minimal "reachable" subset of package and versions that
could appear in pareto-optimal solutions to version resolution for the given set
of required "root" packages, using the following recursive logic:

- A in reqs => A[1] reachable
- A[i] reachable & A[i] depends on B => B[1] reachable
- A[i] reachable & A[i] has reachable conflict => A[i+1] reachable

r is Vec{n}, reachable vector
D is Mat{n,n}, if j depends on P and i = P[1] then on
X is Mat{n,n}, conflicts matrix
S is Mat{n,n}, successor matrix, (i+1, i) with same package

- rule: r′ = r .| (D*r) .| (S*((X*r) .& r))

The function returns a dictionary mapping packages to a vector of PkgVer
structures. Only packages that could be resolved and versions of those packages
that could possibly be chosen by resolution are included.
"""
function reachable(
    pkgs :: Dict{P,PkgInfo{P,V,S}},
    reqs :: Vector{P},
) where {P,V,S}
    # pkg info cache
    pkg_infos = Dict{P,PkgInfo{P,V,S}}()
    pkg!(p) = get!(() -> deps(p), pkg_infos, p)

    # reachable versions (empty)
    reach = Dict{P,Int}()
    queue = [r => 1 for r in reqs]

    # process work queue
    while !empty(queue)
        # get a 
        (p, i) = pop!(queue)
        info = pkg!(p)
        # skip if version doesn't exist
        i ≤ length(info.versions) || continue
        v = info.versions[i]
        
    end
end

function extract(deps::DepsProvider{P,V,S}, reqs::Vector{P}) where {P,V,S}
    # co-compute reachable versions and conflicts between them
    reachable = Pair{P, Int}[p => 1 for p in reqs]
    conflicts = Dict{Int, Set{Int}}()

    # tracking which versions in reachable could conflict
    #   versions: p -> [i] such that
    #       - reachable[i] belongs to package p
    #   interact: p -> [i′] such that
    #       - reachable[i′] depends on p or
    #       - reachable[i′] compat for p
    #
    versions = Dict{P, Vector{Int}}()
    interact = Dict{P, Vector{Int}}()

    # pkg info cache
    cache = Dict{P, PkgInfo{P, V, S}}()
    pkg!(p) = get!(() -> deps(p), cache, p)
    vers!(p) = pkg!(p).vers
    deps!(p) = pkg!(p).deps
    comp!(p) = pkg!(p).comp

    # helper: get versions from reachable
    function get_reachable(i′::Int)
        p′, k′ = reachable[i′]
        v′ = get(vers!(p′), k′, nothing)
        k′, v′, p′
    end

    # helper: check if versions are in reachable
    function in_reachable(p::P, k::Int)
        haskey(versions, p) &&
        any(reachable[j][2] == k for j in versions[p]) ||
        any(reachable[j] == (p => k) for j = i-1:length(reachable))
    end

    # helper: record interactions
    function interact!(p′::P, i::Int)
        interact_p′ = get!(() -> valtype(interact)(), interact, p′)
        i in interact_p′ || push!(interact_p′, i)
    end

    # helper: record conflicts
    function conflict!(i₁, k₁, i₂, k₂; q = i)
        @assert i₁ < i₂
        conflicts_i₁ = get!(() -> Set{Int}(), conflicts, i₁)
        i₂ in conflicts_i₁ && return
        conflicts_i₂ = get!(() -> Set{Int}(), conflicts, i₂)
        # add successor versions to reachable set
        for (i, k) in (i₁ => k₁ + 1, i₂ => k₂ + 1)
            p = reachable[i][1]
            in_reachable(p, k) && continue
            # only add real version or dummy if required package
            if k ≤ length(vers!(p)) + (p ∈ reqs)
                push!(reachable, p => k)
            end
        end
        # record the conflict
        push!(conflicts_i₁, i₂)
        push!(conflicts_i₂, i₁)
    end

    # process enqueued versions
    i = 0
    while i < length(reachable)
        DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

        # take a version off the queue
        k, v, p = get_reachable(i += 1)

        # check if another package has conflict with this one
        haskey(interact, p) && for i′ in interact[p]
            k′, v′, p′ = get_reachable(i′)
            # v′ either depends on p or has compat entry for p
            #   distinguish by whether v === nothing or not
            if v === nothing
                # check if v′ depends on p: can be a false alarm
                # if v′ has compat for p but doesn't depend on p
                haskey(deps!(p′), v′) && p in deps!(p′)[v′] || continue
            else
                c_p′ = comp!(p′)
                # continue if v′ has no compat
                (c_p′v′ = get(c_p′, v′, nothing)) === nothing && continue
                # continue if v′ compat doesn't mention p
                (c_p′v′p = get(c_p′v′, p, nothing)) === nothing && continue
                # continue if v′ compat for p includes v
                v in c_p′v′p && continue
            end
            # conflict:
            #   - v′ depends on p and v === nothing ==> conflict
            #   - v′ compat excluding v !== nothing ==> conflict
            conflict!(i′, k′, i, k)
        end

        # check if this package has conflict with another one
        haskey(comp!(p), v) && for (p′, c_pvp′) in comp!(p)[v]
            interact!(p′, i) # compat is an interaction
            haskey(versions, p′) && for i′ in versions[p′]
                k′, v′ = get_reachable(i′)
                v′ === nothing && continue
                v′ in c_pvp′ && continue
                conflict!(i′, k′, i, k)
            end
        end

        # process dependencies
        haskey(deps!(p), v) && for p′ in deps!(p)[v]
            interact!(p′, i) # dependency is an interaction
            versions_p′ = get(versions, p′, nothing)
            if versions_p′ === nothing
                # new package, insert into reachable
                push!(reachable, p′ => 0, p′ => 1)
                continue
            end
            # dependency => dummy versions are incompatible
            for i′ in (versions_p′[begin], versions_p′[end])
                k′, v′ = get_reachable(i′)
                v′ === nothing || continue
                conflict!(i′, k′, i, k)
            end
        end

        # remember versions of p
        push!(get!(() -> valtype(versions)(), versions, p), i)
    end
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

    # deduplicate versions by their package, dummy flag & conflict set
    kept = collect(reverse(1:length(reachable)))
    while true
        # map all dups to first reachable version
        seen = Dict{Tuple{P, Bool, Set{Int}}, Int}()
        for i in kept
            k, v, p = get_reachable(i)
            x = get!(() -> Set{Int}(), conflicts, i)
            seen[(p, isnothing(v), x)] = i
        end
        length(seen) == length(kept) && break
        keep = Set(values(seen))
        filter!(conflicts) do (i, c)
            filter!(in(keep), c)
            i in keep
        end
        filter!(in(keep), kept)
    end
    reverse!(kept)
    @assert issorted(kept)
    reachable = reachable[kept]
    inds = Dict{Int, Int}(kept .=> 1:length(kept))

    # compute final result structures
    vertices = Pair{P, Union{V, Nothing}}[
        p => get(vers!(p), k, nothing) for (p, k) in reachable
    ]
    conflicts = Set{Int}[
        Set{Int}(inds[j] for j in conflicts[i]) for i in kept
    ]

    # deduplicate final conflict sets
    seen = Dict{Set{Int}, Set{Int}}()
    for (i, x) in enumerate(conflicts)
        conflicts[i] = get!(seen, x, x)
    end

    return vertices, conflicts
end
