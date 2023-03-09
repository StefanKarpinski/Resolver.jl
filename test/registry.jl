using Pkg: depots1
using Pkg.Registry: RegistryInstance, init_package_info!
using Pkg.Types: stdlibs
using Pkg.Versions: VersionSpec

using Resolver
using Resolver: DepsProvider, PkgInfo

const reg_path = joinpath(depots1(), "registries", "General.toml")
const reg_inst = RegistryInstance(expanduser(reg_path))
const reg_dict = Dict(p.name => p for p in values(reg_inst.pkgs))
const excludes = push!(Set(first.(values(stdlibs()))), "julia")

deps = DepsProvider{String, VersionNumber, VersionSpec}() do pkg::String
    info = init_package_info!(reg_dict[pkg])
    vers = sort!(collect(keys(info.version_info)), rev=true)
    deps = Dict(v => String[] for v in vers)
    comp = Dict(v => Dict{String,VersionSpec}() for v in vers)
    # scan versions and populate deps & compat data
    for v in vers
        for (r, d) in info.deps
            v in r && union!(deps[v], keys(d))
        end
        for (r, c) in info.compat
            v in r && mergewith!(intersect, comp[v], c)
        end
    end
    foreach(sort!, values(deps))
    # scrub out excluded deps (stdlibs, julia itself)
    for d in values(deps)
        setdiff!(d, excludes)
    end
    for c in values(comp), x in excludes
        delete!(c, x)
    end
    # deduplicate data structures to save memory
    for i = 1:length(vers)-1, j = i+1:length(vers)
        v, w = vers[i], vers[j]
        deps[v] == deps[w] && (deps[v] = deps[w])
        comp[v] == comp[w] && (comp[v] = comp[w])
    end
    # return resolver PkgInfo data structure
    PkgInfo{String, VersionNumber, VersionSpec}(vers, deps, comp)
end

reqs = sort!(collect(keys(reg_dict)))
filter!(!endswith("_jll"), reqs)
filter!(!in(excludes), reqs)

pkgs = find_packages(deps, reqs)
reach = find_reachable(pkgs, reqs)
filter_reachable!(pkgs, reach)
@time filter_redundant!(pkgs)
