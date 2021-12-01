using Combinatorics
using Random
using Resolver
using Test

function resolve_brute_force(
    packages  :: AbstractVector{<:AbstractVector{<:Integer}},
    conflicts :: AbstractVector{<:Tuple{Integer,Integer}},
)
    N = length(packages)
    M = sum(length, packages)
    P = prod(length, packages)
    P > 0 && log2(P) ≈ sum(log2∘length, packages) ||
        throw(ArgumentError("brute force only works for small problems"))

    # generate all conflict-free solutions
    solution = zeros(Int, N)
    solutions = typeof(solution)[]
    for π in permutations(1:N)
        for S = 0:P-1
            for i = 1:N
                V = length(packages[π[i]])
                S, r = divrem(S, V)
                solution[π[i]] = r
            end
            b = 1
            for i = 1:N
                solution[i] += b
                b += length(packages[i])
            end
            if !any(v₁ ∈ solution && v₂ ∈ solution for (v₁, v₂) in conflicts)
                solution ∉ solutions && push!(solutions, copy(solution))
                break
            end
        end
        isempty(solutions) && break
    end

    # return sorted vector of sorted solutions
    return sort!(solutions)
end

function gen_conflicts(N::I, V::I, C::Integer) where {I<:Integer}
    conflicts = Tuple{I,I}[]
    for p₁ = 0:N-2, p₂ = p₁+1:N-1
        p = N*(N-1)÷2 - (N-p₁)*(N-p₁-1)÷2 + (p₂-p₁) - 1
        @assert 0 ≤ p < N*(N-1)÷2
        X = C >> (p*V^2)
        v = trailing_zeros(X)
        while v < V^2
            v₁, v₂ = divrem(v, V)
            push!(conflicts, (p₁*V + v₁ + 1, p₂*V + v₂ + 1))
            v += 1; v += trailing_zeros(X >> v)
        end
    end
    return conflicts
end
