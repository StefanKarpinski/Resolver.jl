using Test
using Resolver

function solutions_brute_force(
    N::Integer, # number of packages
    V::Integer, # number of versions per package
    conflicts::AbstractVector{<:Tuple{Integer,Integer}},
)
    log(V, V^N) ≈ N ||
        throw(ArgumentError("brute force only works with small problems ($N, $V)"))

    solutions = Vector{Int}[]
    for S = 0:V^N-1
        solution = [V*p + ((S ÷ V^p) % V) + 1 for p = 0:N-1]
        any(v₁ ∈ solution && v₂ ∈ solution for (v₁, v₂) in conflicts) && continue
        push!(solutions, solution)
    end

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

    return solutions
end
