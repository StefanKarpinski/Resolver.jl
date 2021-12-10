module Resolver

export resolve

function resolve(
    packages  :: AbstractVector{<:AbstractVector{<:Integer}},
    conflicts :: AbstractVector{<:Tuple{Integer,Integer}};
    sorted    :: Bool = true,
)
    # counts & sizes
    M = length(packages)                    # number of packages
    N = mapreduce(maximum, max, packages)   # number of versions

    # check packages
    P = zeros(Int, N)
    for (p, versions) in enumerate(packages), v in versions
        1 ≤ v ≤ N ||
            throw(ArgumentError("invalid version index: $v"))
        P[v] == 0 ||
            throw(ArgumentError("version $v in multiple packages"))
        P[v] = p
    end
    for v = 1:N
        P[v] > 0 ||
            throw(ArgumentError("version $v not in any package"))
    end

    # check conflicts
    for (v1, v2) in conflicts
        1 ≤ v1 ≤ N  || throw(ArgumentError("invalid version index: $v1"))
        1 ≤ v2 ≤ N  || throw(ArgumentError("invalid version index: $v2"))
    end

    # no packages, empty solution
    M == 0 && return push!(solutions, Int[])
    # no versions, no solutions
    N == 0 && return solution

    # adjacency lists: incompatible & compatible
    X = [Int[] for v = 1:N]
    C = [Int[] for v = 1:N]
    for (p1, V1) in enumerate(packages), v1 in V1,
        (p2, V2) in enumerate(packages), v2 in V2
        x = p1 == p2 || (v1, v2) ∈ conflicts || (v2, v1) ∈ conflicts
        push!((x ? X : C)[v1], v2)
    end

    # level vector, solution vector, solutions set
    L = ones(Int, N)
    S = zeros(Int, M)
    solutions = typeof(S)[]

    function search!(r::Int = 1, d::Int = 1)
        # find a pivot vertex
        # i = 0
        # for j = 1:N
        #     # check subgraph inclusion at recursion level
        #     L[j] == r || continue
        #     i = j
        #     break
        # end
        # i == 0 && return
        # consider each vertex in pivot set
        for j in 1:N
            # check subgraph inclusion at recursion level
            L[j] == r || continue
            # advance dominance frontier
            d′ = d
            while d′ ≤ length(solutions) && j < solutions[d′][P[j]]
                d′ += 1
            end
            # check advancement of dominance
            d < d′ || length(solutions) < d || continue
            # record version choice
            S[r] = j
            # @show r, L, S[1:r]
            if r == M
                # we have a complete solution
                if d′ > length(solutions)
                    # and it's optimal
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
            L[j] = r-1
        end
    end
    search!()

    # return sorted vector of sorted solutions
    if sorted
        foreach(sort!, solutions)
        sort!(solutions)
    end
    return solutions
end

end # module
