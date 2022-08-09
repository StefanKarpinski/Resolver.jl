using Test
using Random
using Resolver
using Resolver: resolve_core, SetOrVector, DepsProvider, PkgInfo

function resolve_brute_force(
    versions  :: AbstractVector,
    conflicts :: SetOrVector{<:NTuple{2,Integer}};
    optimal   :: Bool = true,
    relax     :: Integer = 0,
)
    # package indices
    P = zeros(Int, length(versions))
    let indices = Dict{eltype(versions), Int}()
        for (v, p) in enumerate(versions)
            P[v] = get!(indices, p, length(indices) + 1)
        end
    end
    perm = sortperm(P)

    # package counts
    M = maximum(P)
    C = zeros(Int, M)
    for p in P
        C[p] += 1
    end

    # number of possible solutions
    ∏ = prod(C)
    log2(∏) ≈ sum(log2, C) ||
        throw(ArgumentError("brute force only works for small problems"))

    # generate all conflict-free solutions
    solution = zeros(Int, M)
    solutions = Vector{Int}[]
    for S = 0:∏-1
        b = 1
        for i = 1:M
            S, r = divrem(S, C[i])
            solution[i] = perm[b + r]
            b += C[i]
        end
        x = sum(v₁ ∈ solution && v₂ ∈ solution for (v₁, v₂) in conflicts; init = 0)
        x ≤ relax && push!(solutions, copy(solution))
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
    foreach(sort!, solutions)
    return sort!(solutions)
end

function gen_conflicts(M::I, V::I, C::Integer) where {I<:Integer}
    conflicts = Tuple{I,I}[]
    for p₁ = 0:M-2, p₂ = p₁+1:M-1
        p = M*(M-1)÷2 - (M-p₁)*(M-p₁-1)÷2 + (p₂-p₁) - 1
        @assert 0 ≤ p < M*(M-1)÷2
        X = C >> (p*V^2)
        v = trailing_zeros(X)
        while v < V^2
            v₁, v₂ = divrem(v, V)
            push!(conflicts, (p₁*V + v₁ + 1, p₂*V + v₂ + 1))
            v += 1; v += trailing_zeros(X >> v)
        end
    end
    return sort!(conflicts)
end

function randu128(k::Integer)
    u = typemax(UInt128)
    while k > 0
        u &= rand(UInt128)
        k -= 1
    end
    return u
end

function make_deps(
    versions  :: Dict{P, <:AbstractVector{V}};
    deps      :: Dict{Pair{P, V}, Vector{P}} =
                 Dict{Pair{P, V}, Vector{P}}(),
    conflicts :: Dict{Tuple{P, P}, Vector{Tuple{V, V}}} =
                 Dict{Tuple{P, P}, Vector{Tuple{V, V}}}(),
) where {P, V}
    A = Dict(p′ => Set(versions[p′]) for p′ in keys(versions))
    DepsProvider{P, V, Set{V}}() do p::P
        vers = versions[p]
        D = Dict(v => d for ((p′, v), d) in deps if p′ == p)
        C = Dict(v => delete!(copy(A), p) for v in vers)
        for ((p₁, p₂), X) in conflicts, (v₁, v₂) in X
            p₁ == p && delete!(C[v₁][p₂], v₂)
            p₂ == p && delete!(C[v₂][p₁], v₁)
        end
        for (v, c) in C
            filter!(c) do (p′, s′)
                s′ ⊈ versions[p′]
            end
        end
        filter!(C) do (v, c)
            !isempty(c)
        end
        PkgInfo{P, V, Set{V}}(collect(vers), D, C)
    end
end
