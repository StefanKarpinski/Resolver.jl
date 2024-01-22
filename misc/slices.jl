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

include("../test/registry.jl")

## load package uuid => name map from registry

function load_packaage_uuid_name_map()
    reg_path = let p = joinpath(depots1(), "registries", "General.toml")
        isfile(p) ? p : splitext(p)[1]
    end
    reg_inst = RegistryInstance(reg_path)
    Dict(string(p.uuid) => p.name for p in values(reg_inst.pkgs))
end

const names = load_packaage_uuid_name_map()

## download package download stats

function load_package_download_stats()
    file = "tmp/package_requests.csv.gz"
    url = "https://julialang-logs.s3.amazonaws.com/public_outputs/current/package_requests.csv.gz"
    isfile(file) || Downloads.download(url, file)
    df = CSV.read(`gzcat $file`, DataFrame)
    filter!(r -> r.status === 200 && isequal(r.client_type, "user"), df)
    @assert allunique(df.package_uuid)
    Dict(names[r.package_uuid] => r.request_addrs for r in eachrow(df))
end

const addrs = load_package_download_stats()
popularity(p) = -get(addrs, p, 0)

## find best installable version of each package

function compute_best_versions(
    sat :: Resolver.SAT{P},
) where {P}
    best_file = "tmp/best_vers.jls"
    if ispath(best_file)
        best = open(deserialize, best_file)
    else
        @info "Computing best pkg versions..."
        best = Dict{P,Int}()
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
    return best
end

## installable packages above median popularity, sorted by popularity

function select_popular_packages(
    best :: Dict{P,Int},
    top  :: Integer,
) where {P}
    packages = sort!(collect(keys(best)))
    sort!(packages, by = popularity)
    N = searchsortedlast(packages, packages[top], by = popularity)
    resize!(packages, N)
end

## pairwise conflicts between best versions

function explicit_conflicts(
    packages :: Vector{P},
    best :: Dict{P,Int},
    sat :: Resolver.SAT{P},
) where {P}
    @info "Constructing explicit conflict graph..."
    N = length(packages)
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
    return G
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

## compute vertex coloring of package graph using RLF algorithm
# https://en.wikipedia.org/wiki/Recursive_largest_first_algorithm

function color_sat_rlf(
    packages :: Vector{P},
    sat :: Resolver.SAT{P},
    G :: AbstractMatrix{Bool},
) where {P}
    @info "Computing coloring (RLF)..."
    N = length(packages)

    # "friendlies" matrix, aka "enemies of enemies", ie two-hop reachability
    H = ((G*G .> 0) .- G .- I(N)) .> 0

    colors = BitVector[]
    U = trues(N) # uncolored nodes

    while any(U)
        # generate a color
        with_temp_clauses(sat) do
            # color data
            C = falses(N)  # current color (independent set)
            n = 0          # color size

            # heuristic data
            A = copy(U)    # available nodes (uncolored & compatible)
            F = fill(0, N) # friendly counts
            X = G*U        # conflict counts (in remaining graph)

            function add_node(i::Int)
                # @assert !S[i]
                # add to color
                C[i] = true
                n += 1
                # add to sat instance
                p = packages[i]
                sat_add(sat, p, best[p])
                sat_add(sat)
                # update heuristic data
                Gᵢ = G[:,i]
                eachnz(Gᵢ) do j, _
                    A[j] = false
                end
                X .-= Gᵢ
                F .+= H[:,i]
                # @assert !any(A .& (G*C .> 0))
                # @assert F == H*C
                # @assert X == G*(U - C)
            end

            # first node maximizes conflicts
            i = max_ind(i->U[i], X, findfirst(U))[1]
            A[i] = false
            add_node(i)

            # progress
            a = count(U)
            prog = Progress(desc="Slice $(length(colors)+1)", a)

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
                i = findfirst(A)
                i === nothing && break
                # maximize friendliness with color
                i, x, tied = max_ind(i->A[i], F, i)
                if tied
                    # minimize external conflicts
                    i = min_ind(X, i) do i
                        A[i] && F[i] == x
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
                A[i] = false
            end

            # progress done
            finish!(prog; showvalues = [
                ("avail", a - n),
                ("count", n),
            ])

            # save color
            push!(colors, C)

            # remove from uncolored nodes
            U .&= .!C

            @assert all(==(1), reduce(+, colors, init=U))
        end
    end

    @info "Colors: $(length(colors))"
    return colors
end

# maximally expand slices, deversifying coverage somewhat
function expand_colors(
    packages :: Vector{P},
    sat :: Resolver.SAT{P},
    colors :: Vector{BitVector},
) where {P}
    @info "Expanding colors into maximal slices..."
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

# compute the implicit conflict graph, i.e. which optimal versions of packages
# are in practice pairwise uninstallable, as opposed to just explicit conflicts
function implicit_conflicts(
    packages :: Vector{P},
    best :: Dict{P,Int},
    slices :: Vector{BitVector},
) where {P}
    @info "Computing implicit conflict graph..."
    N = length(packages)
    # summarize packages by slices they belong to
    S = [[s[i] for s in slices] for i=1:N]
    # for each inclusiion pattern, compute complement of union of slices
    Q = Dict(x => findall(.!reduce(.|, slices[x])) for x in unique(S))
    # compute the true incompatibility graph
    prog = Progress(desc="Implicit conflicts", N)
    G = spzeros(Bool, N, N)
    for (i, p) in enumerate(packages)
        next!(prog; showvalues = [
            ("package", p),
        ])
        for j in Q[S[i]]
            G[i, j] && continue
            q = packages[j]
            sat_assume(sat, p, best[p])
            sat_assume(sat, q, best[q])
            is_satisfiable(sat) && continue
            G[i, j] = true
        end
    end
    return G
end

## top-level code

@info "Loading pkg info..."
info = Resolver.pkg_info(registry.provider())
sat = Resolver.SAT(info)
best = compute_best_versions(sat)

# select top ~1024 most popular packages
top = 1024
packages = select_popular_packages(best, top)
@assert length(packages) ≥ top

# pare down the SAT problem to only popular packages and deps
Resolver.filter_pkg_info!(info, packages)
sat = Resolver.SAT(info)

G₀ = explicit_conflicts(packages, best, sat)
colors₀ = color_sat_rlf(packages, sat, G₀)
slices₀ = expand_colors(packages, sat, colors₀)

G₁ = implicit_conflicts(packages, best, slices₀)
colors₁ = color_sat_rlf(packages, sat, G₁)
slices₁ = expand_colors(packages, sat, colors₁)
