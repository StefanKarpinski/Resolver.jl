module Resolver

const SetOrVector{T} = Union{AbstractSet{T}, AbstractVector{T}}

function resolve(
    compat   :: Function, # (p₁ => v₁, p₂ => v₂) -> Bool
    versions :: AbstractDict{P, <:AbstractVector{V}},
    required :: SetOrVector{P},
) where {P,V}
    vertices, conflicts = vertices_and_conflicts(
        compat,
        versions,
        required,
    )
    packages = P[p for (p, v) in vertices]
    solutions = resolve_core(packages, conflicts)
    # TODO: handle empty solution set
    resolved = Vector{Pair{P,V}}[
        [vertices[v] for v in S if vertices[v][2] !== nothing]
        for S in solutions
    ]
    # only keep the most satisfying solutions
    sat = [sum(p in required for (p, v) in S; init=0) for S in resolved]
    resolved = resolved[sat .≥ maximum(sat)]
end

function vertices_and_conflicts(
    compat   :: Function, # (p₁ => v₁, p₂ => v₂) -> Bool
    versions :: AbstractDict{P, <:AbstractVector{V}},
    required :: SetOrVector{P},
) where {P,V}
    # check compat callback function signature
    hasmethod(compat, Tuple{Pair{P,V}, Pair{P,V}}) &&
    hasmethod(compat, Tuple{Pair{P,V}, Pair{P,Nothing}}) ||
        throw(ArgumentError("compat: callback takes two package-version pairs"))

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

    # co-compute reachable versions and conflicts between them
    package_names = sort!(collect(keys(versions)))
    reachable = [p => Int(p in required) for p in package_names]
    conflicts = Set{NTuple{2,Int}}()
    while true
        clean = true
        for (i₁, (p₁, k₁)) in enumerate(reachable),
            (i₂, (p₂, k₂)) in enumerate(reachable)
            p₁ < p₂ || continue
            v₁ = get(versions[p₁], k₁, nothing)
            v₂ = get(versions[p₂], k₂, nothing)
            (v₁ === nothing || compat(p₁ => v₁, p₂ => v₂)) &&
            (v₂ === nothing || compat(p₂ => v₂, p₁ => v₁)) && continue
            push!(conflicts, minmax(i₁, i₂))
            for (p, k) in (p₁ => k₁ + 1, p₂ => k₂ + 1)
                if k ≤ length(versions[p]) + (p ∈ required) && (p => k) ∉ reachable
                    push!(reachable, p => k)
                    clean = false
                end
            end
        end
        clean && break
    end
    vertices = [p => get(versions[p], k, nothing) for (p, k) in reachable]

    return vertices, conflicts
end

function resolve_core(
    packages  :: AbstractVector,
    conflicts :: SetOrVector{<:NTuple{2,Integer}},
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

    # compatible adjacency lists
    C = [UInt32[] for v = 1:length(packages)]
    for (v₁, p₁) in enumerate(packages),
        (v₂, p₂) in enumerate(packages)
        compatible = p₁ ≠ p₂ && (v₁, v₂) ∉ conflicts && (v₂, v₁) ∉ conflicts
        compatible && push!(C[v₁], v₂)
    end

    # deduplicate nodes by adjacency list
    keep = UInt32[]
    seen = Set{Tuple{eltype(packages), Vector{UInt32}}}()
    for (v, p) in enumerate(packages)
        (p, C[v]) in seen && continue
        push!(seen, (p, C[v]))
        push!(keep, v)
    end
    C = C[keep]
    let kept = Dict(map(reverse, enumerate(keep)))
        for V in C
            filter!(v -> haskey(kept, v), V)
            map!(v -> kept[v], V, V)
        end
    end
    packages = packages[keep]

    # package indices
    P = zeros(UInt32, length(packages))
    let indices = Dict{eltype(packages), UInt32}()
        for (v, p) in enumerate(packages)
            P[v] = get!(indices, p, length(indices) + 1)
        end
    end

    # find all optimal solutions
    solutions = optimal_solutions(P, C)

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

function optimal_solutions(P::Vector{T}, C::Vector{Vector{T}}) where {T<:Integer}
    # number of packages & versions
    M = maximum(P, init=0)
    N = length(C)

    # level vector, solution vector, solutions set
    L = ones(T, N)
    S = zeros(T, M)
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

    foreach(sort!, solutions)
    return sort!(solutions)
end

end # module
