module Slices

export
    popularity,
    compute_best_versions,
    select_popular_packages,
    color_sat_rlf,
    color_sat_dsatur,
    solve_slice,
    expand_slices,
    slices_support,
    conflicts,
    explicit_conflicts,
    implicit_conflicts

using CSV
using DataFrames
using Downloads
using LinearAlgebra
using ProgressMeter
using Serialization
using SparseArrays

using Resolver:
    Resolver,
    resolve_core,
    sat_assume,
    sat_assume_var,
    sat_add,
    sat_add_var,
    is_satisfiable,
    sat_new_variable,
    with_temp_clauses,
    SetOrVec

import Pkg: depots1
import Pkg.Registry: RegistryInstance, init_package_info!

# load package uuid => name map from registry

function load_packaage_uuid_name_map()
    reg_path = let p = joinpath(depots1(), "registries", "General.toml")
        isfile(p) ? p : splitext(p)[1]
    end
    reg_inst = RegistryInstance(reg_path)
    Dict(string(p.uuid) => p.name for p in values(reg_inst.pkgs))
end

const NAMES = load_packaage_uuid_name_map()

# download package download stats

function load_package_download_stats()
    file = joinpath(@__DIR__, "../tmp/package_requests.csv.gz")
    url = "https://julialang-logs.s3.amazonaws.com/public_outputs/current/package_requests.csv.gz"
    isfile(file) || Downloads.download(url, file)
    df = CSV.read(`gzcat $file`, DataFrame)
    filter!(r -> r.status === 200 && isequal(r.client_type, "user"), df)
    @assert allunique(df.package_uuid)
    Dict(NAMES[r.package_uuid] => r.request_addrs for r in eachrow(df))
end

const ADDRS = load_package_download_stats()
popularity(p::AbstractString) = -get(ADDRS, p, 0)
popularity((p, v)::Pair{<:AbstractString,<:Integer}) = popularity(p)

# find best installable version of each package

function compute_best_versions(
    sat :: Resolver.SAT{P},
) where {P}
    best_file = joinpath(@__DIR__, "../tmp/best_vers.jls")
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
    vertices = sort!(collect(best))
    sort!(vertices, by = popularity)
    N = searchsortedlast(vertices, vertices[top], by = popularity)
    resize!(vertices, N)
end

# pairwise conflicts between best versions

function explicit_conflicts(
    vertices :: Vector{Pair{P,Int}},
    info :: Dict{String,<:Resolver.PkgInfo{P}}
) where {P}
    @info "Constructing explicit conflict graph..."
    N = length(vertices)
    G = spzeros(Bool, N, N)
    for (i, (p, v)) in enumerate(vertices)
        T_p = info[p].interacts
        X_p = info[p].conflicts
        for (j, (q, w)) in enumerate(vertices)
            if p == q  && v != w ||
                q in keys(T_p) && X_p[v, T_p[q]+w]
                G[i, j] = true
            end
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
    vertices :: Vector{Pair{P,Int}},
    sat :: Resolver.SAT{P},
    G :: AbstractMatrix{Bool} = conflicts(vertices, sat),
) where {P}
    @info "Computing coloring (RLF)..."
    N = length(vertices)

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
                p, v = vertices[i]
                sat_add(sat, p, v)
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
                ("package", vertices[i][1]),
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
                p, v = vertices[i]
                sat_assume(sat, p, v)
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
    vertices :: Vector{Pair{P,Int}},
    sat :: Resolver.SAT{P},
    G :: AbstractMatrix{Bool} = conflicts(vertices, sat),
) where {P}
    @info "Computing coloring (DSatur)..."
    N = length(vertices)

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
            p, v = vertices[i]
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

function solve_slice(
    vertices :: Vector{Pair{P,Int}},
    sat :: Resolver.SAT{P},
    slice :: BitVector,
) where {P}
    with_temp_clauses(sat) do
        fixed = vertices[slice]
        for (p, v) in fixed
            sat_add(sat, p, v)
            sat_add(sat)
        end
        reqs = unique(first.(fixed))
        only(resolve_core(sat, reqs; max=1, by=popularity))
    end
end

# use resolve to maximally expand slices
function expand_slices!(
    vertices :: Vector{Pair{P,Int}},
    sat :: Resolver.SAT{P},
    slices :: Vector{BitVector},
) where {P}
    @info "Expanding slices..."
    for slice in slices
        sol = solve_slice(vertices, sat, slice)
        for (i, (p, v)) in enumerate(vertices)
            slice[i] = get(sol, p, 0) == v
        end
    end
    return slices
end

function expand_slices(
    vertices :: Vector{Pair{P,Int}},
    sat :: Resolver.SAT{P},
    slices :: Vector{BitVector},
) where {P}
    expand_slices!(vertices, sat, map(copy, slices))
end

function slices_support(
    vertices :: Vector{Pair{P,Int}},
    sat :: Resolver.SAT{P},
    slices :: Vector{BitVector},
) where {P}
    vertices′ = mapreduce(union!, slices) do slice
        sort!(collect(solve_slice(vertices, sat, slice)))
    end
    sort!(vertices′)
    sort!(vertices′, by=popularity)
end

# compute the implicit conflict graph, i.e. which optimal versions of packages
# are in practice pairwise uninstallable, as opposed to just explicit conflicts
function implicit_conflicts(
    vertices :: Vector{Pair{P,Int}},
    sat :: Resolver.SAT{P},
    slices :: Vector{BitVector},
    skip :: SetOrVec{Int} = Set{Int}(),
) where {P}
    @info "Computing implicit conflict graph..."
    N = length(vertices)
    # summarize vertices by slices they belong to
    S = [[s[i] for s in slices] for i=1:N]
    # for each inclusion pattern, compute complement of union of slices
    Q = Dict(x => findall(.!reduce(.|, slices[x])) for x in unique(S))
    # compute the true incompatibility graph
    skip = Set(skip)
    prog = Progress(desc="Implicit conflicts", N)
    G = spzeros(Bool, N, N)
    for (i, (p, v)) in enumerate(vertices)
        next!(prog; showvalues = [
            ("package", p),
        ])
        for j in Q[S[i]]
            if G[i, j] || i ∈ skip && j ∈ skip
                continue
            end
            q, w = vertices[j]
            sat_assume(sat, p, v)
            sat_assume(sat, q, w)
            is_satisfiable(sat) && continue
            G[i, j] = true
        end
    end
    return G
end

function conflicts(
    vertices :: Vector{Pair{P,Int}},
    sat :: Resolver.SAT{P},
) where {P}
    G = explicit_conflicts(vertices, sat.info)
    slices = color_sat_dsatur(vertices, sat, G)
    expand_slices!(vertices, sat, slices)
    implicit_conflicts(vertices, sat, slices)
end

end # module
