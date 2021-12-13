module Resolver

const SetOrVector{T} = Union{AbstractSet{T}, AbstractVector{T}}

function conflicts(
    compatible       :: Function, # ((p1, v1), (p2, v2)) -> Bool
    package_versions :: AbstractDict{<:AbstractString, <:AbstractVector},
)
    # check package versions data structure
    isempty(package_versions) &&
        throw(ArgumentError("packages: no packages"))
    package_names = sort!([String(k) for k in keys(package_versions)])
    for p in package_names
        isempty(package_versions[p]) &&
            throw(ArgumentError("packages: package $p has no versions"))
    end

    # co-compute reachable versions and conflicts between them
    reachable = [(p, 1) for p in package_names]
    conflicts = Set{NTuple{2,Int}}()
    while true
        clean = true
        for (i1, (p1, k1)) in enumerate(reachable),
            (i2, (p2, k2)) in enumerate(reachable)
            p1 < p2 || continue
            v1 = package_versions[p1][k1]
            v2 = package_versions[p2][k2]
            compatible((p1, v1), (p2, v2)) && continue
            push!(conflicts, (i1, i2))
            if k1 < length(package_versions[p1]) && (p1, k1+1) ∉ reachable
                push!(reachable, (p1, k1+1))
                clean = false
            end
            if k2 < length(package_versions[p2]) && (p2, k2+1) ∉ reachable
                push!(reachable, (p2, k2+1))
                clean = false
            end
        end
        clean && break
    end
    versions = [(p, package_versions[p][k]) for (p, k) in reachable]

    return versions, conflicts
end

function resolve_core(
    versions  :: AbstractVector{Package},
    conflicts :: SetOrVector{<:NTuple{2,Integer}},
) where {Package <: Any}
    # check conflicts
    for (v1, v2) in conflicts
        v1 in keys(versions) ||
            throw(ArgumentError("conflicts: invalid version index: $v1"))
        v2 in keys(versions) ||
            throw(ArgumentError("conflicts: invalid version index: $v2"))
        v1 != v2 ||
            throw(ArgumentError("conflicts: package $(versions[v1]): self-conflict $v1"))
        versions[v1] != versions[v2] ||
            throw(ArgumentError("conflicts: package $(versions[v1]): conflict between $v1, $v2"))
    end

    # compatible adjacency lists
    C = [UInt32[] for v = 1:length(versions)]
    for (v1, p1) in enumerate(versions),
        (v2, p2) in enumerate(versions)
        compatible = p1 ≠ p2 && (v1, v2) ∉ conflicts && (v2, v1) ∉ conflicts
        compatible && push!(C[v1], v2)
    end

    # deduplicate nodes by adjacency list
    keep = UInt32[]
    let seen = Set{Tuple{eltype(versions), Vector{UInt32}}}()
        for (v, p) in enumerate(versions)
            (p, C[v]) in seen && continue
            push!(seen, (p, C[v]))
            push!(keep, v)
        end
    end
    C = C[keep]
    let kept = Dict(map(reverse, enumerate(keep)))
        for V in C
            filter!(v -> haskey(kept, v), V)
            map!(v -> kept[v], V, V)
        end
    end
    versions = versions[keep]

    # package indices
    P = zeros(UInt32, length(versions))
    let indices = Dict{eltype(versions), UInt32}()
        for (v, p) in enumerate(versions)
            P[v] = get!(indices, p, length(indices) + 1)
        end
    end

    # find all optimal solutions
    solutions = find_optimal_solutions(P, C)

    # translate back to original version indices
    Vector{Int}[Int[keep[v] for v in S] for S in solutions]
end

function find_optimal_solutions(P::Vector{T}, C::Vector{Vector{T}}) where {T<:Integer}
    # number of packages & versions
    M = maximum(P)
    N = length(C)

    # level vector, solution vector, solutions set
    L = ones(T, N)
    S = zeros(T, M)
    solutions = typeof(S)[]

    function search!(r::Int = 1, d::Int = 1)
        for j in 1:N
            # check subgraph inclusion at this recursion level
            L[j] == r || continue
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
            # restrict next subgraph to neighbors
            for k in C[j]
                L[k] += (L[k] == r)
            end
            # recursion
            search!(r + 1, d′)
            # restore the levels vector
            for k in C[j]
                L[k] = min(L[k], r)
            end
            # exclude vertex from future iterations
            L[j] = r - 1
        end
    end
    search!()

    foreach(sort!, solutions)
    return sort!(solutions)
end

end # module
