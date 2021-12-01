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
