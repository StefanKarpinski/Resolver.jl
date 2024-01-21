using CSV
using DataFrames
using Downloads
using LinearAlgebra
using ProgressMeter
using Resolver
using Serialization
using SparseArrays

using Resolver:
    sat_assume,
    sat_add,
    is_satisfiable,
    extract_solution!,
    optimize_solution!,
    with_temp_clauses

import Pkg: depots1
import Pkg.Registry: RegistryInstance, init_package_info!

## load package uuid => name map from registry

reg_path = let p = joinpath(depots1(), "registries", "General.toml")
    isfile(p) ? p : splitext(p)[1]
end
reg_inst = RegistryInstance(reg_path)
names = Dict(string(p.uuid) => p.name for p in values(reg_inst.pkgs))

## download package download stats

file = "tmp/package_requests.csv.gz"
url = "https://julialang-logs.s3.amazonaws.com/public_outputs/current/package_requests.csv.gz"
isfile(file) || Downloads.download(url, file)
df = CSV.read(`gzcat $file`, DataFrame)
filter!(r -> r.status === 200 && isequal(r.client_type, "user"), df)
@assert allunique(df.package_uuid)
const addrs = Dict(names[r.package_uuid] => r.request_addrs for r in eachrow(df))
popularity(p) = -get(addrs, p, 0)

## load resolution problem

include("../test/registry.jl")
rp = registry.provider()
info = Resolver.pkg_info(rp)

## construct SAT problem

sat = Resolver.SAT(info)
@assert is_satisfiable(sat)
@assert !is_satisfiable(sat, keys(info))

# type parameters
(P, V) = (String, VersionNumber)

## find best installable version of each package
best_file = "tmp/best_vers.jls"
if ispath(best_file)
    const best = open(deserialize, best_file)
else
    const best = Dict{P,Int}()
    let prog = Progress(desc="Best versions", length(sat.info))
        for p in keys(sat.info)
            next!(prog; showvalues = [
                ("package", p),
            ])
            for i = 1:length(sat.info[p].versions)
                sat_assume(sat, p, i)
                if is_satisfiable(sat)
                    best[p] = i
                    break
                end
            end
        end
    end
    open(best_file, write=true) do io
        serialize(io, best)
    end
end

## installable packages above median popularity, sorted by popularity

const packages = sort!(collect(keys(best)))
sort!(packages, by = popularity)
N = searchsortedlast(packages, packages[1024], by = popularity)
resize!(packages, N)

# pare down the SAT problem as well

Resolver.filter_pkg_info!(info, packages)
sat = Resolver.SAT(info)

## pairwise conflicts between best versions

G = spzeros(Bool, N, N)
for (i, p) in enumerate(packages)
    v = best[p]
    info_p = sat.info[p]
    for (q, b) in info_p.interacts
        w = get(best, q, 0)
        w == 0 && continue
        info_p.conflicts[v, b+w] || continue
        j = findfirst(==(q), packages)
        j === nothing && continue
        G[i, j] = true
    end
end

## index search helpers

# max index search
# - checks for ties
# - no early return
function max_ind(
    p  :: Function,       # index predicate
    v  :: AbstractVector, # values vector
    i₀ :: Int,            # initial index
)
    @assert p(i₀)
    max_i = i₀
    max_x = v[i₀]
    tied = false
    for i = i₀+1:length(v)
        p(i) || continue
        x = v[i]
        if x > max_x
            max_i = i
            max_x = x
        elseif x == max_x
            tied = true
        end
    end
    return max_i, max_x, tied
end

# min index search
# - doesn't check for ties
# - returns early if min_x hits zero
function min_ind(
    p  :: Function,       # index predicate
    v  :: AbstractVector, # values vector
    i₀ :: Int,            # initial index
)
    @assert p(i₀)
    min_i = i₀
    min_x = v[i₀]
    for i = i₀+1:length(v)
        min_x ≤ 0 && break
        p(i) || continue
        x = v[i]
        if x < min_x
            min_i = i
            min_x = x
        end
    end
    return min_i
end

# helper for iterating non-zeros of a sparse vector
function eachnz(f::Function, S::SparseVector)
    for (i, j) in enumerate(S.nzind)
        v = S.nzval[i]
        !iszero(v) && f(j, v)
    end
end

## compute package slices

function color_sat_rlf(
    G :: AbstractMatrix{Bool},
    packages :: Vector{P} = packages,
    sat :: Resolver.SAT{P} = sat,
)
    N = length(packages)

    # "friendlies" matrix, aka "enemies of enemies", ie two-hop reachability
    H = ((G*G .> 0) .- G .- I(N)) .> 0

    slices = BitVector[]
    A = trues(N) # available nodes

    while any(A)
        # generate a slice
        with_temp_clauses(sat) do
            # slice data
            S = falses(N)  # current slice (independent set)
            n = 0          # slice size

            # heuristic data
            C = copy(A)    # candidate nodes (compatible & available)
            F = fill(0, N) # friendly counts
            X = G*A        # conflict counts (in remaining graph)

            function add_node(i::Int)
                # @assert !S[i]
                # add to slice
                S[i] = true
                n += 1
                # add to sat instance
                p = packages[i]
                sat_add(sat, p, best[p])
                sat_add(sat)
                # update heuristic data
                Gᵢ = G[:,i]
                eachnz(Gᵢ) do j, _
                    C[j] = false
                end
                X .-= Gᵢ
                F .+= H[:,i]
                # @assert !any(C .& (G*S .> 0))
                # @assert F == H*S
                # @assert X == G*(A - S)
            end

            # first node maximizes conflicts
            i = max_ind(i->A[i], X, findfirst(A))[1]
            C[i] = false
            add_node(i)

            # progress
            a = count(A)
            prog = Progress(desc="Slice $(length(slices)+1)", a)

            # progress update
            next!(prog; showvalues = [
                ("avail", a - n),
                ("count", n),
                ("index", i),
                ("package", packages[i]),
                ("sat", true),
            ])

            # add rest of nodes
            while true
                i = findfirst(C)
                i === nothing && break
                # maximize friendliness with slice
                i, x, tied = max_ind(i->C[i], F, i)
                if tied
                    # minimize external conflicts
                    i = min_ind(X, i) do i
                        C[i] && F[i] == x
                    end
                end
                # check for satisfiability (not pure graph coloring)
                p = packages[i]
                sat_assume(sat, p, best[p])
                s = is_satisfiable(sat)
                s && add_node(i)
                # progress update
                next!(prog; showvalues = [
                    ("avail", a - n),
                    ("count", n),
                    ("index", i),
                    ("package", p),
                    ("sat", s),
                ])
                # don't consider again
                C[i] = false
            end

            # progress done
            finish!(prog; showvalues = [
                ("avail", a - n),
                ("count", n),
            ])

            # save slice
            push!(slices, S)

            # remove from available nodes
            A .&= .!S

            @assert all(==(1), reduce(+, slices, init=A))
        end
    end

    return slices
end

colors = color_sat_rlf(G)

# There may be versions without explicit incompatibilities that cannot actually
# be used together. The naive way of computing this is to cheack for pair that
# isn't explicitly incompatible if the SAT problem is satisfiable. This way too
# slow, however, as it's O(N^2) and each SAT call is non-trivial.
#
# One way of cutting that down is to use the explicit incompatibility graph and
# graph coloring to partition the graph into disjoint compatible subsets. Since
# this uses the explicit incompatibility graph for the heuristics, it may not
# produce the optimal results, but we know each slice is compatible and any true
# incompatibilities must be between the slices. Next, we grow those slices as
# large as we can. If S is the original disjoint slice, write S* for the maximal
# set "grown from" S. The pairs we must consider are the pairs that appear in
# none of the slices. When expanding slices, we'll want to preferrentially add
# versions that we haven't already added to other slices. A good heuristic might
# be to add packages that belong to the fewest expanded slices already.
#
# Consider p ∈ S. We're looking for q that are not in any slice with p. If q is
# in some slice with p, then we know they're compatible. The values of q that we
# don't need to consider are all of those that appear in some slice with p, i.e.
# the union of all the slices that p appears in. The complement of that union is
# the set of q that we need to consider for incompatibility. Write Q(p) for the
# complement of the union of expanded slices that contain p. We want to consider
# pairs of the form p, q ∈ Q(p) when looking for potentially incompatible pairs.
#
# Once we've done this to compute the true pairwise compatibility graph, we can
# run our slice coloring algorithm again and perhaps get better results since
# our input graph for heuristics is more accurate.

# maximally expand slices, deversifying coverage somewhat
function expand_colors(
    colors :: Vector{BitVector},
    packages :: Vector{P} = packages,
    sat :: Resolver.SAT{P} = sat,
)
    slices = map(copy, colors)
    order = copy(packages)
    todos = Set(packages)
    for slice in slices
        with_temp_clauses(sat) do
            for (i, p) in enumerate(packages)
                slice[i] || continue
                sat_add(sat, p, best[p])
                sat_add(sat)
            end
            for p in order
                i = findfirst(==(p), packages)
                slice[i] && continue
                sat_assume(sat, p, best[p])
                is_satisfiable(sat) || continue
                sat_add(sat, p, best[p])
                sat_add(sat)
                slice[i] = true
                delete!(todos, p)
            end
        end
        sort!(order, by = !in(todos))
    end
    return slices
end

slices = expand_colors(colors)

# summarize packages by slices they belong to
sinc = [[s[i] for s in slices] for i=1:N]
rinc = Dict{Vector{Bool},Vector{Int}}()
for (i, k) in enumerate(sinc)
    push!(get!(()->Int[], rinc, k), i)
end

# for each pattern, compute the complement of union of its slices
Q = Dict(x => findall(.!reduce(.|, slices[x])) for x in keys(rinc))

# compute the true incompatibility graph
G′ = copy(G)
for (i, p) in enumerate(packages)
    for j in Q[sinc[i]]
        G′[i, j] && continue
        q = packages[j]
        sat_assume(sat, p, best[p])
        sat_assume(sat, q, best[q])
        is_satisfiable(sat) && continue
        G′[i, j] = true
    end
end

colors′ = color_sat_rlf(G′)
slices′ = expand_colors(colors′)
