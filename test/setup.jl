using Test
using Resolver

function resolve_brute_force(
    packages  :: AbstractVector{<:AbstractVector{<:Integer}},
    conflicts :: AbstractVector{<:Tuple{Integer,Integer}},
)
    N = length(packages)
    M = sum(length, packages)
    Π = prod(length, packages)
    Π > 0 && log2(Π) ≈ sum(log2∘length, packages) ||
        throw(ArgumentError("brute force only works for small problems"))

    # generate all conflict-free solutions
    solution = zeros(Int, N)
    solutions = Vector{Int}[]
    for S = 0:Π-1
        b = 1
        for i = 1:N
            V = length(packages[i])
            S, r = divrem(S, V)
            solution[i] = b + r
            b += V
        end
        any(v₁ ∈ solution && v₂ ∈ solution for (v₁, v₂) in conflicts) && continue
        push!(solutions, copy(solution))
    end

    # filter out sub-optimal solutions
    filter!(solutions) do s₁
        all(solutions) do s₂
            equal = true
            for (v₁, v₂) in zip(s₁, s₂)
                v₁ < v₂ && return true
                equal &= (v₁ == v₂)
            end
            return equal
        end
    end

    # return sorted vector of sorted solutions
    return sort!(map(sort!, solutions))
end
