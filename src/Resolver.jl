module Resolver

export resolve

function resolve(
    packages  :: AbstractVector{<:AbstractVector{<:Integer}},
    conflicts :: AbstractVector{<:Tuple{Integer,Integer}},
)
    # counts & sizes
    M = length(packages)                    # number of packages
    N = mapreduce(maximum, max, packages)   # number of versions

    # check packages
    P = zeros(Int, N)
    for (p, V) in enumerate(packages)
        length(V) > 0 ||
            throw(ArgumentError("packages: package $p has no versions"))
        for v in V
            1 ≤ v ≤ N ||
                throw(ArgumentError("packages: invalid version index: $v"))
            P[v] == 0 ||
                throw(ArgumentError("packages: version $v in multiple packages"))
            P[v] = p
        end
    end
    # check that we haven't skipped any version number
    for v = 1:N
        P[v] > 0 ||
            throw(ArgumentError("packages: version $v not in any package"))
    end

    # check conflicts
    for (v1, v2) in conflicts
        1 ≤ v1 ≤ N ||
            throw(ArgumentError("conflicts: invalid version index: $v1"))
        1 ≤ v2 ≤ N ||
            throw(ArgumentError("conflicts: invalid version index: $v2"))
    end

    # no packages, empty solution
    M > 0 || return [Int[]]


    # compatible adjacency lists
    C = [UInt32[] for v = 1:N]
    for (p1, V1) in enumerate(packages), v1 in V1,
        (p2, V2) in enumerate(packages), v2 in V2
        if p1 ≠ p2 && (v1, v2) ∉ conflicts && (v2, v1) ∉ conflicts
            push!(C[v1], v2)
        end
    end

    # deduplicate nodes by adjacency list
    keep = UInt32[]
    let seen = Set{Tuple{UInt32,Vector{UInt32}}}()
        for (p, V) in enumerate(packages), v in V
            (p, C[v]) in seen && continue
            push!(seen, (p, C[v]))
            push!(keep, v)
        end
    end
    P = P[keep]
    C = C[keep]
    let d = Dict(map(reverse, enumerate(keep)))
        for V in C
            filter!(V) do v
                haskey(d, v)
            end
            map!(V, V) do v
                d[v]
            end
        end
    end
    N = length(keep)

    # level vector, solution vector, solutions set
    L = ones(UInt32, N)
    S = zeros(UInt32, M)
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
    sort!(solutions)
    Vector{Int}[Int[keep[v] for v in S] for S in solutions]
end

end # module
