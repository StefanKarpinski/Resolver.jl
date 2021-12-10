using Test
using Random
using Resolver

function resolve_brute_force(
    packages  :: AbstractVector{<:AbstractVector{<:Integer}},
    conflicts :: AbstractVector{<:Tuple{Integer,Integer}};
    optimal   :: Bool = true,
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
    optimal && filter!(solutions) do s₁
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

function randu128(k::Integer)
    u = typemax(UInt128)
    while k > 0
        u &= rand(UInt128)
        k -= 1
    end
    return u
end
