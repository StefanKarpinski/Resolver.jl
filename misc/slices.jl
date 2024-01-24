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
    sat_assume_var,
    sat_add,
    sat_add_var,
    is_satisfiable,
    sat_new_variable,
    extract_solution!,
    solution,
    optimize_solution!,
    with_temp_clauses

import Pkg: depots1
import Pkg.Registry: RegistryInstance, init_package_info!

include("../test/registry.jl")

# load package uuid => name map from registry

function load_packaage_uuid_name_map()
    reg_path = let p = joinpath(depots1(), "registries", "General.toml")
        isfile(p) ? p : splitext(p)[1]
    end
    reg_inst = RegistryInstance(reg_path)
    Dict(string(p.uuid) => p.name for p in values(reg_inst.pkgs))
end

const names = load_packaage_uuid_name_map()

# download package download stats

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

# find best installable version of each package

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

# installable packages above median popularity, sorted by popularity

function select_popular_packages(
    best :: Dict{P,Int},
    top  :: Integer,
) where {P}
    packages = sort!(collect(keys(best)))
    sort!(packages, by = popularity)
    N = searchsortedlast(packages, packages[top], by = popularity)
    resize!(packages, N)
end

# pairwise conflicts between best versions

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

# optimal (min/max) index search

# NOTE: passing a stop_val invalidates `tied`
# in cases where the early stop_val is found

for (name, lt) in ((:min_ind, :(<)), (:max_ind, :(>)))
    @eval function $name(
        predicate :: Function,       # index predicate
        values    :: AbstractVector, # values vector
        start_ind :: Int;            # start index
        stop_val  :: Any = nothing,  # stop value
    )
        @assert predicate(start_ind)
        opt_i = start_ind
        opt_x = values[start_ind]
        tied = false
        for i = start_ind+1:length(values)
            opt_x == stop_val && break
            predicate(i) || continue
            x = values[i]
            if $lt(x, opt_x)
                opt_i = i
                opt_x = x
            elseif x == opt_x
                tied = true
            end
        end
        return opt_i, opt_x, tied
    end
    @eval $name(
        predicate :: AbstractVector{Bool},
        values    :: AbstractVector,
        start_ind :: Int = findfirst(predicate);
        stop_val  :: Any = nothing,
    ) = $name(i->predicate[i], values, start_ind; stop_val)
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
            i = max_ind(U, X)[1]
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
                i, x, tied = max_ind(A, F, i)
                if tied
                    # minimize external conflicts
                    i, = min_ind(X, i, stop_val=0) do i
                        A[i] && F[i] == x
                    end
                end
                # check for satisfiability (not done for pure graph coloring)
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

## compute vertex coloring of package graph using DSatur algorithm
# https://en.wikipedia.org/wiki/Recursive_largest_first_algorithm

function color_sat_dsatur(
    packages :: Vector{P},
    sat :: Resolver.SAT{P},
    G :: AbstractMatrix{Bool},
) where {P}
    @info "Computing coloring (DSatur)..."
    N = length(packages)

    colors = BitVector[]
    vars = Int[] # SAT variables for colors

    U = trues(N)                   # uncolored nodes
    A = SparseVector{Bool,Int64}[] # adjacency vectors
    S = fill(0, N)                 # saturation counts
    X = G*U                        # conflicts counts

    with_temp_clauses(sat) do
        while (i = findfirst(U)) !== nothing
            i, s, tied = max_ind(U, S, i)
            if tied
                # maximize external conflicts
                i, = max_ind(X, i) do i
                    U[i] && S[i] == s
                end
            end
            p = packages[i]
            v = best[p]
            Gᵢ = G[:,i]
            # add node to first color it's compatible with
            success = false
            for j in sortperm(colors, by = count)
                c = vars[j]
                # check satisfiability: c & p@v
                sat_assume_var(sat, c)
                sat_assume(sat, p, v)
                if is_satisfiable(sat)
                    # add to existing color
                    colors[j][i] = true
                    # add to color: c => p@v
                    sat_add_var(sat, -c)
                    sat_add(sat, p, v)
                    sat_add(sat)
                    # update heuristic data
                    A[j] .|= Gᵢ
                    eachnz(Gᵢ) do k, _
                        S[k] = sum(a[k] for a in A)
                    end
                    # end loop successfully
                    success = true
                    break
                end
            end
            if !success
                # create new color
                push!(colors, falses(N))
                colors[end][i] = true
                # new SAT variable
                c = sat_new_variable(sat)
                push!(vars, c)
                sat_add_var(sat, -c)
                sat_add(sat, p, v)
                sat_add(sat)
                # new adjacency vector
                push!(A, Gᵢ)
                S .+= Gᵢ
            end
            # remove node
            X .-= Gᵢ
            U[i] = false
            # @assert X == G*U
            # @assert S == [sum(a[k] for a in A) for k=1:N]
            # @assert all(==(1), reduce(+, colors, init=U))
        end
    end

    @info "Colors: $(length(colors))"
    return colors
end

# maximally expand slices, deversifying coverage somewhat
function expand_slices(
    packages :: Vector{P},
    sat :: Resolver.SAT{P},
    slices :: Vector{BitVector},
) where {P}
    @info "Expanding slices..."
    slices = map(copy, slices)
    order = copy(packages)
    todos = Set(packages)
    N = length(packages)
    prog = Progress(desc="Expanding slices", N*length(slices))
    for (k, slice) in enumerate(slices)
        n = length(slice)
        if n < N
            resize!(slice, N)
            slice[n+1:N] .= false
        end
        with_temp_clauses(sat) do
            for (i, p) in enumerate(packages)
                slice[i] || continue
                sat_add(sat, p, best[p])
                sat_add(sat)
                next!(prog, showvalues = [
                    ("slice", k),
                    ("package", p),
                ])
            end
            for p in order
                i = findfirst(==(p), packages)
                slice[i] && continue
                next!(prog, showvalues = [
                    ("slice", k),
                    ("package", p),
                ])
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

# sat add a slice

function slice_dict(
    packages :: Vector{P},
    best :: Dict{P,Int},
    slice :: AbstractVector{Bool},
) where {P}
    vers = Dict{P,Int}()
    for (i, p) in enumerate(packages)
        i ≤ length(slice) || break
        slice[i] || continue
        vers[p] = best[p]
    end
    return vers
end

## top-level code

@info "Loading pkg info..."
info = Resolver.pkg_info(registry.provider())
sat = Resolver.SAT(info)
best = compute_best_versions(sat)

# select top most popular packages
top = 1024
packages = select_popular_packages(best, top)
@assert length(packages) ≥ top

# pare down the SAT problem to only popular packages and deps
Resolver.filter_pkg_info!(info, packages)
sat = Resolver.SAT(info)

G₀ = explicit_conflicts(packages, best, sat)
colors₀ = color_sat_dsatur(packages, sat, G₀)
slices₀ = expand_slices(packages, sat, colors₀)

G = implicit_conflicts(packages, best, slices₀)
colors = color_sat_dsatur(packages, sat, G)
slices = expand_slices(packages, sat, colors)

# check slices
for slice in slices
    vers = slice_dict(packages, best, slice)
    sat_assume(sat, vers)
    @assert is_satisfiable(sat)
end

# turn slices into solutions
sols = [slice_dict(packages, best, slice) for slice in slices]
for sol in sols
    with_temp_clauses(sat) do
        for (p, v) in sol
            sat_add(sat, p, v)
            sat_add(sat)
        end
        for p in packages
            p in keys(sol) && continue
            sat_assume(sat, p)
            is_satisfiable(sat) || continue

            # require some version of p
            sat_add(sat, p)
            sat_add(sat)
            extract_solution!(sat, sol)

            # optimize p's version
            best[p] < sol[p] &&
            optimize_solution!(sat, sol) do
                for i = 1:sol[p]-1
                    sat_add(sat, p, i)
                end
                sat_add(sat)
            end

            # fix p's version
            sat_add(sat, p, sol[p])
            sat_add(sat)
        end
    end
end

vers = Dict{String,Vector{Int}}()
for sol in sols
    for (p, v) in sol
        V = get!(()->Int[], vers, p)
        v in V || push!(V, v)
    end
end
map(sort!, values(vers))
