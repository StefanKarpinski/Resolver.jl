# Pkg's module structure is bonkers, wrap in a single module
module Pkg
    import Pkg: depots1
    import Pkg.Registry: RegistryInstance, init_package_info!
    import Pkg.Types: stdlibs
    import Pkg.Versions: VersionSpec
end
import .Pkg

const REG_PATH = joinpath(Pkg.depots1(), "registries", "General.toml")
const EXCLUDES = push!(Set(first.(values(Pkg.stdlibs()))), "julia")

function registry_provider(
    reg_path :: AbstractString = REG_PATH;
    excludes :: SetOrVec{<:AbstractString} = EXCLUDES,
)
    reg_inst = Pkg.RegistryInstance(reg_path)
    reg_dict = Dict(p.name => p for p in values(reg_inst.pkgs) if p.name ∉ excludes)

    DepsProvider{String, VersionNumber, Pkg.VersionSpec}(keys(reg_dict)) do pkg::String
        info = Pkg.init_package_info!(reg_dict[pkg])
        vers = sort!(collect(keys(info.version_info)), rev=true)
        deps = Dict(v => String[] for v in vers)
        comp = Dict(v => Dict{String,Pkg.VersionSpec}() for v in vers)
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
        # deduplicate data structures to save memory
        for i = 1:length(vers)-1, j = i+1:length(vers)
            v, w = vers[i], vers[j]
            deps[v] == deps[w] && (deps[v] = deps[w])
            comp[v] == comp[w] && (comp[v] = comp[w])
        end
        # return resolver PkgData structure
        PkgData{String, VersionNumber, Pkg.VersionSpec}(vers, deps, comp)
    end
end
