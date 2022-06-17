module Resolver

DEBUG::Bool = false

const SetOrVector{T} = Union{AbstractSet{T}, AbstractVector{T}}

include("DepProviders.jl")

function resolve(deps::DepProvider{P,V}, reqs::SetOrVector{P}) where {P,V}
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
            sat = [sum(p in required for (p, v) in S; init=0) for S in resolved]
            max_sat = maximum(sat; init = 1-isempty(vertices))
            resolved = resolved[sat .≥ max_sat]
            !isempty(resolved) && break
        end
        resolved
    end
end

function prepare(deps::DepCache{P,V,S}, reqs::SetOrVector{P}) where {P,V,S}
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

    # co-compute reachable versions and conflicts between them
    reachable = Pair{P, Int}[p => 1 for p in reqs]
    conflicts = Set{Tuple{Int, Int}}()

    # tracking which requirements could conflict
    interact = Dict{P, Dict{P, Vector{Int}}}()
    versions = Dict{P, Vector{Int}}()

    # accessors: versions, dependencies & compat
    vers(p) = versions!(deps, p) :: Vector{V}
    deps(p) = depends!(deps, p)  :: Dict{V, Vector{P}}
    comp(p) = compat!(deps, p)   :: Dict{V, Dict{P, S}}

    # helper to record conflicts
    function conflict!(i₁, k₁, i₂, k₂)
        x = minmax(i₁, i₂)
        x in conflicts && return
        push!(conflicts, x)
        for (p, k) in (p₁ => k₁ + 1, p₂ => k₂ + 1)
            # continue if we already have p => k
            any(reachable[j][2] == k for j in versions[p]) && continue
            # only add real version or dummy for required package
            if k ≤ length(vers(p)) + (p ∈ reqs)
                push!(reachable, p => k)
            end
        end
    end

    # process enqueued versions
    i = 0
    while i < length(reachable)
        DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

        # take a version off the queue
        p, k = reachable[i += 1]
        v = vers(p)[k]

        # make sure dependencies are reachable
        for d in deps(p)[v]
            haskey(versions, d) && continue
            push!(reachable, d => 0, d => 1)
        end

        # check if another package has conflict with this one
        for (p′, I′) in interact[p]
            c_p′ = comp(p′)
            for i′ in I′
                p_, k′ = reachable[i′]
                @assert p′ == p_
                v′ = vers(p′)[k′]
                (c_p′v′ = get(c_p′, v′, nothing)) === nothing && continue
                (c_p′v′p = get(c_p′v′, p, nothing)) === nothing && continue
                v in c_p′v′p && continue
                conflict!(i, k, i′, k′)
            end
        end
        # check if this package has conflict with another one
        for (p′, c_pvp′) in comp(p)[v]
            for i′ in versions[p′]
                p_, k′ = reachable[i′]
                @assert p′ == p_
                v′ = vers(p′)[k′]
                v′ in c_pvp′ && continue
                conflict!(i, k, i′, k′)
            end
            # remember potential conflicts
            push!(get!(() -> Int[], interact[p′], p), i)
        end
        # ensure interaction dict for p exists
        get!(interact, p) do
            valtype(interact)()
        end
        # remember new version of p
        push!(get!(() -> valtype(versions), versions, p), i)
    end
    DEBUG && println(@__FILE__, ":", @__LINE__, " @ ", time()-t₀)

    # convert version indices back to versions
    vertices = Pair{P,Union{V,Nothing}}[
        p => get(vers(p), k, nothing) for (p, k) in reachable
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
