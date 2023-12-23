module registry

import Resolver: DepsProvider, PkgData

import Pkg: depots1
import Pkg.Registry: RegistryInstance, init_package_info!
import Pkg.Types: stdlibs
import Pkg.Versions: VersionSpec

const REG_PATH = joinpath(depots1(), "registries", "General.toml")
const EXCLUDES = push!(Set(first.(values(stdlibs()))), "julia")

function provider(
    reg_path :: AbstractString = REG_PATH;
    excludes :: AbstractSet{<:AbstractString} = EXCLUDES,
)
    reg_inst = RegistryInstance(reg_path)
    reg_dict = Dict(p.name => p
        for p in values(reg_inst.pkgs) if p.name ∉ excludes)

    DepsProvider(keys(reg_dict)) do pkg :: AbstractString
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
                v in r && mergewith!(∩, comp[v], c)
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
        # deduplicate data structures to save some memory
        for i = 1:length(vers)-1, j = i+1:length(vers)
            v, w = vers[i], vers[j]
            deps[v] == deps[w] && (deps[v] = deps[w])
            comp[v] == comp[w] && (comp[v] = comp[w])
        end
        # return resolver PkgData structure
        PkgData(vers, deps, comp)
    end
end

end # module
