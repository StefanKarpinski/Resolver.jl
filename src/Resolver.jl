module Resolver

export resolve

# for debugging
dd(A::AbstractArray) = map(reverse∘bitstring, permutedims(A))

function resolve(
    packages  :: AbstractVector{<:AbstractVector{<:Integer}},
    conflicts :: AbstractVector{<:Tuple{Integer,Integer}};
    sort      :: Bool = true,
)
    # counts & sizes
    M = length(packages)                    # number of packages
    N = mapreduce(maximum, max, packages)   # number of versions

    # check packages
    let counts = zeros(Int, N)
        for (p, versions) in enumerate(packages), v in versions
            1 ≤ v ≤ N ||
                throw(ArgumentError("invalid version index: $v"))
            counts[v] += 1
        end
        for v = 1:N
            counts[v] == 1 ||
                throw(ArgumentError("version $v in $(counts[v]) packages"))
        end
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
        x = p1 == p2 || (v1, v2) ∈ conflicts || (v2, v2) ∈ conflicts
        push!((x ? X : C)[v1], v2)
    end

    # level vector, solution vector, solutions set
    L = ones(Int, N)
    S = zeros(Int, M)
    solutions = typeof(S)[]

    function search!(r::Int = 1)
        # find a pivot vertex
        i = 0
        for j = 1:N
            l = L[j]
            if l == r
                i = j
                break
            end
        end
        i == 0 && return
        # consider each vertex in pivot set
        for j in X[i]
            l = L[j]
            if l == r
                S[r] = j
                if r == M
                    push!(solutions, copy(S))
                    continue
                end
                for k in C[j]
                    L[k] += (L[k] == r)
                end
                search!(r + 1)
                for k in C[j]
                    L[k] = min(L[k], r)
                end
                L[j] = r-1
            end
        end
    end
    search!()

    # return sorted vector of sorted solutions
    if sort
        foreach(sort!, solutions)
        sort!(solutions)
    end
    return solutions
end

end # module
