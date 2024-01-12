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
best = Dict{P,Int}()
let prog = Progress(length(sat.info), desc="Best versions")
    for p in keys(sat.info)
        for i = 1:length(sat.info[p].versions)
            sat_assume(sat, p, i)
            if is_satisfiable(sat)
                best[p] = i
                break
            end
        end
        next!(prog; showvalues = [
            ("package", p),
        ])
    end
end

# compute package slices

const packages = sort!(collect(keys(best)), by = popularity)
const todo = Set(packages)
const sols = Dict{P,Int}[]
const sol = Dict{P,Int}()

while !isempty(todo)
    # generate the next slice
    slice = length(sols) + 1
    prog = Progress(desc="Slice $slice", 2*length(packages))
    with_temp_clauses(sat) do
        # satisfy as many todos as possible
        for p in sort!(packages, by = !in(todo))
            # update progress
            next!(prog; showvalues = [
                ("package", p),
                ("todos", length(todo)),
            ])

            # check if p@i is feasible
            i = best[p]
            sat_assume(sat, p, i)
            is_satisfiable(sat) || continue

            # require p@i
            sat_add(sat, p, i)
            sat_add(sat)

            # delete from todos
            delete!(todo, p)
        end
        # optimize remaining versions
        for p in sort!(packages, by = popularity)
            # update progress
            next!(prog; showvalues = [
                ("package", p),
                ("todos", length(todo)),
                ("solution", length(sol)),
            ])

            # check if p is feasible
            sat_assume(sat, p)
            is_satisfiable(sat) || continue

            # require some version of p
            sat_add(sat, p)
            sat_add(sat)
            extract_solution!(sat, sol)

            # optimize p's version
            best[p] < sol[p] &&
            optimize_solution!(sat, sol) do
                # check if some strictly better version exists
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
    push!(sols, copy(sol))
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

nodes = Pair{String,Int}[]
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
    sol = Dict{String,Int}()
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
