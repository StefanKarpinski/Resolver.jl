using CSV
using DataFrames
using Downloads
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

include("test/registry.jl")
rp = registry.provider()
info = Resolver.pkg_info(rp)

# order packages by download address count

packages = sort!(collect(keys(info)))
sort!(packages, by = p -> get(addrs, p, 0), rev = true)

# solve SAT problems

(P, V) = (String, VersionNumber)

using Resolver:
    sat_assume,
    sat_add,
    is_satisfiable,
    extract_solution!,
    optimize_solution!

sat = Resolver.SAT(info)
@assert Resolver.is_satisfiable(sat)
@assert !Resolver.is_satisfiable(sat, keys(info))

sol = Dict{P,Int}()
for (k, p) in enumerate(packages)
    sat_assume(sat, p)
    is_satisfiable(sat) || continue

    println("$k $(length(sol)) $p")

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
end
