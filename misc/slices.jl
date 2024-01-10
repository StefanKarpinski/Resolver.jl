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
addrs = Dict(names[r.package_uuid] => r.request_addrs for r in eachrow(df))

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

packages = sort!(collect(keys(best)))
sort!(packages, by = p -> -get(addrs, p, 0))

todo = Set(packages[1:1024])
sols = Dict{P,Int}[]
sol = Dict{P,Int}()
while true
    slice = length(sols) + 1
    # compute an optimized solution
    with_temp_clauses(sat) do
        prog = Progress(desc="Slice $slice", length(packages))
        for (k, p) in enumerate(packages)
            # check if p is feasible
            sat_assume(sat, p)
            is_satisfiable(sat) || continue

            # require some version of p
            sat_add(sat, p)
            sat_add(sat)
            extract_solution!(sat, sol)

            # optimize p's version
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

            sol[p] ≤ best[p] && delete!(todo, p)

            # update progress
            update!(prog, k; showvalues = [
                ("package", p),
                ("solution", length(sol)),
                ("todos", length(todo)),
            ])
        end
    end
    push!(sols, copy(sol))

    # check if we're done
    isempty(todo) && break

    # prioritize not-yet-optimized packages
    sort!(packages, by = p -> -get(addrs, p, 0))
    sort!(packages, by = !in(todo))
end

#=
nodes = Dict{Pair{String,Int},Int}()
let j = 0
    for sol in sols, (p, i) in sol
        haskey(nodes, p => i) && continue
        nodes[p => i] = j += 1
    end
end

using Graphs
G = SimpleGraph(length(nodes))
for sol in sols, (p, i) in sol, (q, j) in sol
    p < q || continue
    add_edge!(G, nodes[p => i], nodes[q => j])
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
