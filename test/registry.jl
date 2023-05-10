using Pkg: depots1
using Pkg.Registry: RegistryInstance, init_package_info!
using Pkg.Types: stdlibs
using Pkg.Versions: VersionSpec

const reg_path = joinpath(depots1(), "registries", "General.toml")
const reg_inst = RegistryInstance(expanduser(reg_path))
const reg_dict = Dict(p.name => p for p in values(reg_inst.pkgs))
const excludes = push!(Set(first.(values(stdlibs()))), "julia")

dp = DepsProvider{String, VersionNumber, VersionSpec}() do pkg::String
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

#=
all_names = sort!(collect(keys(reg_dict)))
filter!(!endswith("_jll"), all_names)
filter!(!in(excludes), all_names)
all_pkgs = find_packages(dp, all_names)
filter_reachable!(all_pkgs, all_names)
filter_redundant!(all_pkgs)
all_ix = find_interacts(all_pkgs)

const pairs = Tuple{String,String}[]
for p in all_names, q in get(all_ix, p, String[])
    p < q || continue
    reqs = [p, q]
    pkgs = find_packages(dp, reqs)
    filter_reachable!(pkgs, reqs)
    filter_redundant!(pkgs)
    all(length(pkgs[p].versions) > 1 for p in reqs) || continue
    push!(pairs, (p, q))
    @show reqs
end
=#

solver = expanduser("~/dev/kissat/build/kissat")

# reqs = String.(split("JSON", ','))
reqs = String.(split("YaoLang,ZXCalculus", ','))
pkgs = find_packages(dp, reqs)
filter_reachable!(pkgs, reqs)
filter_redundant!(pkgs)

problem = gen_sat("tmp/problem.cnf", pkgs, reqs)
run(ignorestatus(`$solver -q $problem`))
