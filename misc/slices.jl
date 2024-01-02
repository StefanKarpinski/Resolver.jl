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

url = "https://julialang-logs.s3.amazonaws.com/public_outputs/current/package_requests.csv.gz"
file = Downloads.download(url)
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
let prog = Progress(length(sat.info), desc="Best")
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
    # compute an optimized solution
    with_temp_clauses(sat) do
        prog = Progress(length(packages),
            desc = "Slice $(length(sols)+1) ($(length(todo)) todos)")
        for (k, p) in enumerate(packages)
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

            # update progress
            update!(prog, k; showvalues = [
                ("package", p),
                ("solution", length(sol)),
            ])
        end
    end
    push!(sols, copy(sol))

    # delete fully optimized packages from todo set
    for (p, i) in sol
        i ≤ best[p] && delete!(todo, p)
    end
    isempty(todo) && break

    # prioritize not-yet-optimized packages
    sort!(packages, by = !in(todo))
end
