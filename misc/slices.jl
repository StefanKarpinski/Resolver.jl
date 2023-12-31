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

using Printf

sat = Resolver.SAT(info)
@assert Resolver.is_satisfiable(sat)
@assert !Resolver.is_satisfiable(sat, keys(info))

reqs = Set{String}()
for (i, p) in enumerate(packages)
    p in reqs && continue
    @printf("%i %0.4f %s\n", length(reqs), i/length(packages), p)
    push!(reqs, p)
    pkgs, vers = resolve(sat, reqs; max=1)
    if isempty(vers)
        delete!(reqs, p)
        continue
    end
    vers = vers[:,1]
    Resolver.fix_versions(sat, pkgs, vers)
    for (p, v) in zip(pkgs, vers)
        v !== nothing && push!(reqs, p)
    end
end
