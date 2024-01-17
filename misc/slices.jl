using CSV
using DataFrames
using Downloads
using ProgressMeter
using Resolver

# load package uuid => name map from registry

import Pkg: depots1
import Pkg.Registry: RegistryInstance, init_package_info!

reg_path = let p = joinpath(depots1(), "registries", "General.toml")
    isfile(p) ? p : splitext(p)[1]
end
reg_inst = RegistryInstance(reg_path)
names = Dict(string(p.uuid) => p.name for p in values(reg_inst.pkgs))

# download package download stats

file = "tmp/package_requests.csv.gz"
url = "https://julialang-logs.s3.amazonaws.com/public_outputs/current/package_requests.csv.gz"
isfile(file) || Downloads.download(url, file)
df = CSV.read(`gzcat $file`, DataFrame)
filter!(r -> r.status === 200 && isequal(r.client_type, "user"), df)
@assert allunique(df.package_uuid)
const addrs = Dict(names[r.package_uuid] => r.request_addrs for r in eachrow(df))
popularity(p) = -get(addrs, p, 0)

# load resolution problem

include("../test/registry.jl")
rp = registry.provider()
info = Resolver.pkg_info(rp)

# construct SAT problem

using Resolver:
    sat_assume,
    sat_add,
    is_satisfiable,
    extract_solution!,
    optimize_solution!,
    with_temp_clauses

sat = Resolver.SAT(info)
@assert is_satisfiable(sat)
@assert !is_satisfiable(sat, keys(info))

# type parameters
(P, V) = (String, VersionNumber)

# find best installable version of each package
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

# installable packages, sorted by popularity

const packages = sort!(collect(keys(best)))
sort!(packages, by = popularity)

# pairwise compatibility of best versions

const N = length(best)
const X = falses(N, N)

for (k, p) in enumerate(packages)
    i = best[p]
    info_p = sat.info[p]
    for (q, b) in info_p.interacts
        j = get(best, q, 0)
        j == 0 && continue
        info_p.conflicts[i, b+j] || continue
        l = findfirst(==(q), packages)
        X[k, l] = true
    end
end

# compute package slices

const sats = typeof(sat)[]
const slices = Dict{P,Int}[]
const conflicts = Vector{Int}[]

prog = Progress(desc="Slices", N)
for i = 1:N # (1, 5093, 2, 3)
    p = packages[i]
    # update progress
    showvalues = [("package", p), ("index", i)]
    for s = 1:length(slices)
        push!(showvalues, ("slice $s", length(slices[s])))
    end
    next!(prog; showvalues)
    # actual slice computation
    best_p = best[p]
    candidates = Int[]
    for (s, sat) in enumerate(sats)
        sat_assume(sat, p, best_p)
        is_satisfiable(sat) && push!(candidates, s)
    end
    if isempty(candidates) # add to new slice
        push!(sats, Resolver.SAT(info))
        push!(slices, Dict{P,Int}())
        push!(conflicts, fill(0, N))
        push!(candidates, length(sats))
    end
    @assert length(candidates) > 0
    if length(candidates) == 1
        s = only(candidates)
    else # multiple candidate slices
        # pick the slice that leaves as many optimal versions of the
        # remaining packages compatible with some existing slice
        best_s = best_c = -1
        best_z = maximum(length, slices) + 1
        for s in candidates
            conflicts[s] .+= X[:, i]
            c = sum(any(C[j] == 0 for C in conflicts) for j=i+1:N)
            conflicts[s] .-= X[:, i]
            z = length(slices[s])
            if c > best_c || c == best_c && z < best_z
                best_s = s
                best_c = c
                best_z = z
            end
        end
        s = best_s
    end
    # add to chosen slice
    sat_add(sats[s], p, best_p)
    sat_add(sats[s])
    slices[s][p] = best_p
    conflicts[s] .+= X[:, i]
end

#=
# # only use versions from sols after this
# for p in packages
#     vers = Int[]
#     for sol in sols
#         i = get(sol, p, nothing)
#         i isa Int || continue
#         i in vers && continue
#         push!(vers, i)
#     end
#     sort!(vers)
#     for i in vers
#         sat_add(sat, p, i)
#     end
#     sat_add(sat)
# end

if !@isdefined(common)
    common = Dict(reduce(∩, sols))
end
for sol in sols, p in keys(common)
    delete!(sol, p)
end

nodes = Pair{P,Int}[]
for sol in sols, (p, i) in sol
    (p => i) in nodes && continue
    push!(nodes, p => i)
end
rnodes = Dict(map(reverse, enumerate(nodes)))

n = length(nodes)
G = BitMatrix(undef, n, n)
fill!(G, false)
for sol in sols, (p, i) in sol, (q, j) in sol
    G[rnodes[p => i], rnodes[q => j]] = true
end

using GraphModularDecomposition
T = sort!(StrongModuleTree(G))

# fill in all possible edges
let prog = Progress(desc="Pairs", size(G,1)),
    sol = Dict{P,Int}()
    for k = 1:size(G,1)
        p, i = nodes[k]
        with_temp_clauses(sat) do
            sat_add(sat, p, i)
            sat_add(sat)
            for l = k+1:size(G,2)
                G[k, l] && continue
                q, j = nodes[l]
                sat_assume(sat, q, j)
                is_satisfiable(sat) || continue
                extract_solution!(sat, sol)
                for (p′, i′) in sol
                    k′ = rnodes[p′ => i′]
                    for (q′, i′) in sol
                        p′ < q′ || continue
                        l′ = rnodes[q′ => j′]
                        G[k′, l′] = G[l′, k′] = true
                    end
                end
            end
        end
    end
end

=#

#=
using Graphs
G = SimpleGraph(length(nodes))
for sol in sols, (p, i) in sol, (q, j) in sol
    p < q || continue
    add_edge!(G, rnodes[p => i], rnodes[q => j])
end

function is_cograph(G::Graph)
    if length(vertices(G)) <= 1
        return true
    end
    if !is_connected(G)
        for V in connected_components(G)
            is_cograph(complement(G[V])) || return false
        end
        return true
    end
    Ḡ = complement(G)
    if !is_connected(Ḡ)
        for V in connected_components(Ḡ)
            is_cograph(G[V]) || return false
        end
        return true
    end
    return false
end

is_cograph(G)
=#
