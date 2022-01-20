module Resolver

const SetOrVector{T} = Union{AbstractSet{T}, AbstractVector{T}}
const DepsType = Union{
    Function, # Pair{P,V} --> SetOrVector{P}
    AbstractDict{Pair{P,V}, <:SetOrVector{P}} where {P,V},
}

const NO_CONFLICTS = (::Pair, ::Pair) -> true
const NO_DEPS = ((::Pair{P},) where P) -> Set{P}()

function resolve(
    versions :: AbstractDict{P, <:AbstractVector{V}},
    required :: SetOrVector{P};
    compat   :: Function = NO_CONFLICTS, # (p₁ => v₁, p₂ => v₂) -> Bool
    deps     :: DepsType = NO_DEPS,      # (p₁ => v₁) -> SetOrVector{P}
) where {P,V}
    if deps isa AbstractDict
        deps_dict = deps
        deps = let no_deps = Set{P}()
            pv::Pair{P,V} -> get(deps_dict, pv, no_deps)
        end
    end
    vertices, conflicts = vertices_and_conflicts(
        versions,
        required;
        compat,
        deps,
    )
    packages = P[p for (p, v) in vertices]
    dummies = Bool[isnothing(v) for (p, v) in vertices]
    let resolved
        M = length(versions)
        for relax = 0:M*(M-1)÷2
            solutions = resolve_core(packages, conflicts; dummies, relax)
            # sort solutions
            for S in solutions
                sort!(S, by = v -> packages[v])
            end
            sort!(solutions)
            # map back to version identifiers
            resolved = Vector{Pair{P,V}}[
                [vertices[v] for v in S if !dummies[v]]
                    for S in solutions
            ]
            # only keep the most satisfying solutions
            sat = [sum(p in required for (p, v) in S; init=0) for S in resolved]
            resolved = resolved[sat .≥ maximum(sat; init = 1-isempty(versions))]
            !isempty(resolved) && break
        end
        resolved
    end
end

function vertices_and_conflicts(
    versions :: AbstractDict{P, <:AbstractVector{V}},
    required :: SetOrVector{P};
    compat   :: Function, # (p₁ => v₁, p₂ => v₂) -> Bool
    deps     :: Function, # (p₁ => v₁) -> SetOrVector{P}
) where {P,V}
    # check compat callback function signature
    hasmethod(compat, Tuple{Pair{P,V}, Pair{P,V}}) ||
        throw(ArgumentError("compat: callback takes two package-version pairs"))

    # check deps callback function signature
    hasmethod(deps, Tuple{Pair{P,V}}) ||
        throw(ArgumentError("deps: callback takes a package-version pair"))

    # check package versions data structure
    let seen = Set{V}()
        for (p, vers) in versions
            isempty(vers) &&
                throw(ArgumentError("versions: package $p: no versions"))
            for v in vers
                v in seen &&
                    throw(ArgumentError("versions: package $p: duplicate version $v"))
                push!(seen, v)
            end
            empty!(seen)
        end
    end

    # check required packages set
    for p in required
        p in keys(versions) ||
            throw(ArgumentError("required: package $p: not in versions"))
    end

    # cache of deps lists
    deps_cache = Dict{Pair{P,V}, Set{P}}()
    deps!(pv::Pair{P,V}) = get!(deps_cache, pv) do
        s = deps(pv)
        s isa AbstractSet ? s : Set(s)
    end

    # co-compute reachable versions and conflicts between them
    package_names = sort!(collect(keys(versions)))
    reachable = Pair{P,Int}[p => Int(p in required) for p in package_names]
    conflicts = Set{NTuple{2,Int}}()
    conflicted = Set(required)
    while true
        clean = true
        for (i₁, (p₁, k₁)) in enumerate(reachable),
            (i₂, (p₂, k₂)) in enumerate(reachable)
            p₁ < p₂ || continue
            v₁ = get(versions[p₁], k₁, nothing)
            v₂ = get(versions[p₂], k₂, nothing)
            if v₁ !== nothing && v₂ !== nothing
                compat(p₁ => v₁, p₂ => v₂) &&
                compat(p₂ => v₂, p₁ => v₁) && continue
            elseif v₁ === v₂ === nothing
                continue
            elseif v₁ === nothing
                p₁ ∈ deps!(p₂ => v₂) || continue
            elseif v₂ === nothing
                p₂ ∈ deps!(p₁ => v₁) || continue
            end
            push!(conflicts, minmax(i₁, i₂))
            for (p, k) in (p₁ => k₁ + 1, p₂ => k₂ + 1)
                if k ≤ length(versions[p]) + (p ∈ required) && (p => k) ∉ reachable
                    push!(reachable, p => k)
                    clean = false
                end
                push!(conflicted, p)
            end
        end
        clean && break
    end
    vertices = Pair{P,Union{V,Nothing}}[
        p => get(versions[p], k, nothing)
            for (p, k) in reachable if p in conflicted
    ]

    return vertices, conflicts
end

function resolve_core(
    packages  :: AbstractVector,
    conflicts :: SetOrVector{<:NTuple{2,Integer}};
    dummies   :: Vector{Bool} = fill(false, length(packages)),
    relax     :: Integer = 0,
)
    # check conflicts
    for (v₁, v₂) in conflicts
        v₁ in keys(packages) ||
            throw(ArgumentError("conflicts: invalid version index: $v₁"))
        v₂ in keys(packages) ||
            throw(ArgumentError("conflicts: invalid version index: $v₂"))
        v₁ != v₂ ||
            throw(ArgumentError("conflicts: package $(packages[v₁]): self-conflict $v₁"))
        packages[v₁] != packages[v₂] ||
            throw(ArgumentError("conflicts: package $(packages[v₁]): conflict between $v₁, $v₂"))
    end

    # conflict count matrix
    N = length(packages)
    X = zeros(UInt32, N, N)
    for (v₁, p₁) in enumerate(packages),
        (v₂, p₂) in enumerate(packages)
        x = (v₁, v₂) ∈ conflicts || (v₂, v₁) ∈ conflicts
        d = x && (dummies[v₁] || dummies[v₂])
        X[v₁, v₂] = X[v₂, v₁] = p₁ == p₂ ? relax + 2 : d ? relax + 1 : x
    end

    # deduplicate nodes by conflict sets
    keep = UInt32[]
    seen = Set{Vector{UInt32}}()
    for (v, p) in enumerate(packages)
        X[:, v] in seen && continue
        push!(seen, X[:, v])
        push!(keep, v)
    end
    X = X[keep, keep]
    packages = packages[keep]

    # package indices
    P = zeros(UInt32, length(packages))
    let indices = Dict{eltype(packages), UInt32}()
        for (v, p) in enumerate(packages)
            P[v] = get!(indices, p, length(indices) + 1)
        end
    end

    # find all optimal solutions
    solutions = optimal_solutions(P, X, Int(relax))

    # translate back to original version indices
    Vector{Int}[Int[keep[v] for v in S] for S in solutions]
end

# The core solver is based on naive k-clique listing algorithm in Chiba &
# Nishizeki 1985, "Arboricity and subgraph listing algorithms". It is also
# similar to the well-known Bron-Kerbosch maximal clique listing algorithm, but
# doesn't require keeping track of the set of excluded nodes because the clique
# size is known in advance. This version lists only optimal solutions in the
# following sense: we have a total ordering on the versions of each package and
# we must choose some version of each package (possibly the special "no version"
# version); this induces a partial ordering on solutions and the solutions that
# are listed are those that are optimal by this partial ordering. These are the
# conflict-free solutions that a user might want to choose based on version
# preferences. If there is a conflict between the latest versions of packages A
# and B, they might want the latest version of A and an older version of B or
# the latest version of B and an older version of A. Whereas the total number of
# feasible conflict-free solutions is exponential in the number of conflict
# edges, the number of optimal solutions tends to be fairly small.
#
# The pivoting optimization that can be used to speed up Bron-Kerbosch cannot be
# used here because if we both limit vertices considered at each recursion level
# to the neighborhood of some node and require choosing a node that dominates a
# specific solution in the set of already-found solutions, then we miss valid,
# optimal solutions. However, the domination requirement acts as a form of
# pivoting anyway: after finding an optimal solution, we only consider versions
# that improve on that solution for some package. Those versions must also
# conflict with that solution since otherwise it wouldn't be optimal. So rather
# than pivoting around a node, we are pivoting around an entire solution.

function optimal_solutions(P::Vector{T}, X::Matrix{T}, c::Int = 0) where {T<:Integer}
    # number of packages & versions
    M = maximum(P, init=0)
    L = maximum(X, init=1) + M
    N = size(X, 1)

    # conflict count vector, level vector, solution vector, solutions set
    C = zeros(T, N)
    S = zeros(T, M)
    solutions = typeof(S)[]
    M == 0 && return push!(solutions, S)

    function search!(c::Int, r::Int = 1, d::Int = 1)
        for j = 1:N
            # check subgraph inclusion at this level
            c′ = c - C[j]
            c′ ≥ 0 || continue
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
            # restrict next subgraph to neighbors without too many conflicts
            for k = 1:N
                C[k] += X[k, j]
            end
            # recursion
            search!(c′, r + 1, d′)
            # restore the conflict count & level vectors
            for k = 1:N
                C[k] -= X[k, j]
            end
            # exclude vertex from future iterations
            C[j] += r*L
        end
        for j = 1:N
            dr = divrem(C[j], L)
            if dr[1] ≥ r
                C[j] = dr[2]
            end
        end
    end
    search!(c)

    foreach(sort!, solutions)
    return sort!(solutions)
end

end # module
