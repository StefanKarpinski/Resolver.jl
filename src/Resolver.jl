module Resolver

DEBUG::Bool = false

const SetOrVector{T} = Union{AbstractSet{T}, AbstractVector{T}}

# Type Parameters:
#  P = package type
#  V = version type
#  S = version set type

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

function resolve(deps::DepsProvider{P,V}, reqs::Vector{P}) where {P,V}
    @eval t₀::Float64 = time()
    vertices, conflicts = prepare(deps, reqs)
    packages = P[p for (p, v) in vertices]
    dummies = Bool[isnothing(v) for (p, v) in vertices]
    let resolved
        M = length(vertices)
        for relax = 0:M*(M-1)÷2
            solutions = resolve_core(packages, conflicts; dummies, relax)
            # sort solutions
            for S in solutions
                sort!(S, by = v -> packages[v])
            end
            sort!(solutions)
            # map back to version identifiers
            resolved = Vector{Pair{P,V}}[
                [vertices[v] for v in S if !dummies[v]]
                    for S in solutions
            ]
            # only keep the most satisfying solutions
            sat = [sum(p in reqs for (p, v) in S; init=0) for S in resolved]
            max_sat = maximum(sat; init = 1-isempty(vertices))
            resolved = resolved[sat .≥ max_sat]
            !isempty(resolved) && break
        end
        resolved
    end
end

function prepare(deps::DepsProvider{P,V,S}, reqs::Vector{P}) where {P,V,S}
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

    # co-compute reachable versions and conflicts between them
    reachable = Pair{P, Int}[p => 1 for p in reqs]
    conflicts = Set{Tuple{Int, Int}}()

    # tracking which requirements could conflict
    interact = Dict{P, Dict{P, Vector{Int}}}()
    versions = Dict{P, Vector{Int}}()

    # pkg info cache
    cache = Dict{P, PkgInfo{P, V, S}}()
    pkg!(p) = get!(() -> deps(p), cache, p)
    vers!(p) = pkg!(p).versions
    deps!(p) = pkg!(p).depends
    comp!(p) = pkg!(p).compat

    function get_reachable(p::P, i::Int)
        p′, k = reachable[i]
        @assert p′ == p
        k, get(vers!(p), k, nothing)
    end

    function in_reachable(p::P, k::Int)
        haskey(versions, p) &&
        any(reachable[j][2] == k for j in versions[p]) ||
        any(reachable[j] == (p => k) for j = i-1:length(reachable))
    end

    # helper to record conflicts
    function conflict!(i₁, k₁, i₂, k₂; q = i)
        @assert i₁ < i₂
        (i₁, i₂) in conflicts && return
        for (i, k) in (i₁ => k₁ + 1, i₂ => k₂ + 1)
            p = reachable[i][1]
            in_reachable(p, k) && continue
            # only add real version or dummy if required package
            if k ≤ length(vers!(p)) + (p ∈ reqs)
                push!(reachable, p => k)
            end
        end
        push!(conflicts, (i₁, i₂))
    end

    # process enqueued versions
    i = 0
    while i < length(reachable)
        DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

        # take a version off the queue
        p, k = reachable[i += 1]
        v = get(vers!(p), k, nothing)

        # make sure dependencies are reachable
        haskey(deps!(p), v) && for p′ in deps!(p)[v]
            haskey(versions, p′) && continue
            push!(reachable, p′ => 0, p′ => 1)
        end

        # check if another package has conflict with this one
        haskey(interact, p) && for (p′, I′) in interact[p]
            c_p′ = comp!(p′)
            for i′ in I′
                k′, v′ = get_reachable(p′, i′)
                # continue if v′ has no compat
                (c_p′v′ = get(c_p′, v′, nothing)) === nothing && continue
                # continue if v′ compat doesn't mention p
                (c_p′v′p = get(c_p′v′, p, nothing)) === nothing && continue
                # continue if v′ compat for p includes v
                v in c_p′v′p && continue
                # conflict: v′ has compat for p and doesn't include v
                conflict!(i′, k′, i, k)
            end
        end

        # check if this package has conflict with another one
        haskey(comp!(p), v) && for (p′, c_pvp′) in comp!(p)[v]
            haskey(versions, p′) && for i′ in versions[p′]
                k′, v′ = get_reachable(p′, i′)
                v′ in c_pvp′ && continue
                conflict!(i′, k′, i, k)
            end
            # remember as a potential conflict
            interact_p′ = get!(() -> valtype(interact)(), interact, p′)
            push!(get!(() -> Int[], interact_p′, p), i)
        end

        # remember new version of p
        push!(get!(() -> valtype(versions)(), versions, p), i)
    end
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

    # convert version indices back to versions
    vertices = Pair{P,Union{V,Nothing}}[
        p => get(vers!(p), k, nothing) for (p, k) in reachable
    ]

    return vertices, conflicts
end

function resolve_core(
    packages  :: AbstractVector,
    conflicts :: SetOrVector{<:NTuple{2,Integer}};
    dummies   :: Vector{Bool} = fill(false, length(packages)),
    relax     :: Integer = 0,
)
    # check conflicts
    for (v₁, v₂) in conflicts
        v₁ in keys(packages) ||
            throw(ArgumentError("conflicts: invalid version index: $v₁"))
        v₂ in keys(packages) ||
            throw(ArgumentError("conflicts: invalid version index: $v₂"))
        v₁ != v₂ ||
            throw(ArgumentError("conflicts: package $(packages[v₁]): self-conflict $v₁"))
        packages[v₁] != packages[v₂] ||
            throw(ArgumentError("conflicts: package $(packages[v₁]): conflict between $v₁, $v₂"))
    end

    # conflict count matrix
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)
    N = length(packages)
    X = zeros(UInt32, N, N)
    for (v₁, p₁) in enumerate(packages),
        (v₂, p₂) in enumerate(packages)
        x = (v₁, v₂) ∈ conflicts || (v₂, v₁) ∈ conflicts
        d = x && (dummies[v₁] || dummies[v₂])
        X[v₁, v₂] = X[v₂, v₁] = p₁ == p₂ ? relax + 2 : d ? relax + 1 : x
    end
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

    # deduplicate nodes by conflict sets
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)
    keep = UInt32[]
    seen = Set{Vector{UInt32}}()
    for (v, p) in enumerate(packages)
        X[:, v] in seen && continue
        push!(seen, X[:, v])
        push!(keep, v)
    end
    X = X[keep, keep]
    packages = packages[keep]
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

    # package indices
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)
    P = zeros(UInt32, length(packages))
    let indices = Dict{eltype(packages), UInt32}()
        for (v, p) in enumerate(packages)
            P[v] = get!(indices, p, length(indices) + 1)
        end
    end
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

    # find all optimal solutions
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)
    solutions = optimal_solutions(P, X, Int(relax))
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

    # translate back to original version indices
    Vector{Int}[Int[keep[v] for v in S] for S in solutions]
end

# The core solver is based on naive k-clique listing algorithm in Chiba &
# Nishizeki 1985, "Arboricity and subgraph listing algorithms". It is also
# similar to the well-known Bron-Kerbosch maximal clique listing algorithm, but
# doesn't require keeping track of the set of excluded nodes because the clique
# size is known in advance. This version lists only optimal solutions in the
# following sense: we have a total ordering on the versions of each package and
# we must choose some version of each package (possibly the special "no version"
# version); this induces a partial ordering on solutions and the solutions that
# are listed are those that are optimal by this partial ordering. These are the
# conflict-free solutions that a user might want to choose based on version
# preferences. If there is a conflict between the latest versions of packages A
# and B, they might want the latest version of A and an older version of B or
# the latest version of B and an older version of A. Whereas the total number of
# feasible conflict-free solutions is exponential in the number of conflict
# edges, the number of optimal solutions tends to be fairly small.
#
# The pivoting optimization that can be used to speed up Bron-Kerbosch cannot be
# used here because if we both limit vertices considered at each recursion level
# to the neighborhood of some node and require choosing a node that dominates a
# specific solution in the set of already-found solutions, then we miss valid,
# optimal solutions. However, the domination requirement acts as a form of
# pivoting anyway: after finding an optimal solution, we only consider versions
# that improve on that solution for some package. Those versions must also
# conflict with that solution since otherwise it wouldn't be optimal. So rather
# than pivoting around a node, we are pivoting around an entire solution.

function optimal_solutions(P::Vector{T}, X::Matrix{T}, c::Int = 0) where {T<:Integer}
    # number of packages & versions
    M = maximum(P, init=0)
    L = maximum(X, init=1) + M
    N = size(X, 1)

    # conflict count vector, level vector, solution vector, solutions set
    C = zeros(T, N)
    S = zeros(T, M)
    solutions = typeof(S)[]
    M == 0 && return push!(solutions, S)

    function search!(c::Int, r::Int = 1, d::Int = 1)
        for j = 1:N
            # check subgraph inclusion at this level
            c′ = c - C[j]
            c′ ≥ 0 || continue
            # record version choice
            S[r] = j
            # advance dominance frontier (if necessary & possible)
            l = length(solutions)
            d′ = d
            if d ≤ l
                # check for advancement
                j < solutions[d][P[j]] || continue
                # advance past already-dominated solutions
                while true
                    d′ += 1
                    d′ ≤ l || break
                    S′ = solutions[d′]
                    any(k < S′[P[k]] for r′ = 1:r for k = S[r′]) || break
                end
            end
            # check for a complete solution
            if r == M
                if d′ > l # it's optimal, save it!
                    push!(solutions, sort(S, by = k -> P[k]))
                end
                break
            end
            # restrict next subgraph to neighbors without too many conflicts
            for k = 1:N
                C[k] += X[k, j]
            end
            # recursion
            search!(c′, r + 1, d′)
            # restore the conflict count & level vectors
            for k = 1:N
                C[k] -= X[k, j]
            end
            # exclude vertex from future iterations
            C[j] += r*L
        end
        for j = 1:N
            dr = divrem(C[j], L)
            if dr[1] ≥ r
                C[j] = dr[2]
            end
        end
    end
    search!(c)

    foreach(sort!, solutions)
    return sort!(solutions)
end

end # module
