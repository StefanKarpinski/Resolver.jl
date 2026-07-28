module registry

import Resolver: DepsProvider, PkgData

import Pkg: depots1
import Pkg.Registry: RegistryInstance, init_package_info!
import Pkg.Types: stdlibs
import Pkg.Versions: VersionSpec

const EXCLUDES = push!(Set(first.(values(stdlibs()))), "julia")
const REG_PATH =
    let p = joinpath(depots1(), "registries", "General.toml")
        isfile(p) ? p : splitext(p)[1]
    end

# Julia 1.14 made the parsed registry data uuid-keyed instead of name-keyed:
#   deps   value: Dict{String,UUID} (<=1.13)  ->  Set{UUID}            (>=1.14)
#   compat value: Dict{String,VersionSpec}     ->  Dict{UUID,VersionSpec}
# This provider is name-based, so recover names via the registry. Stdlib/julia
# uuids aren't registered, so they're skipped -- they'd be excluded anyway.
_dep_names(d::AbstractDict, reg) = keys(d)
_dep_names(d::AbstractSet, reg) =
    (reg.pkgs[u].name for u in d if haskey(reg.pkgs, u))
_compat_names(c::AbstractDict{<:AbstractString}, reg) = c
_compat_names(c::AbstractDict{Base.UUID}, reg) =
    Dict(reg.pkgs[u].name => s for (u, s) in c if haskey(reg.pkgs, u))

# share one object among equal-valued entries; there are only a few
# distinct values per package, so scanning them beats comparing all pairs
function dedup_values!(d::AbstractDict, vers::Vector)
    distinct = valtype(d)[]
    for v in vers
        x = d[v]
        shared = false
        for y in distinct
            if y == x
                d[v] = y
                shared = true
                break
            end
        end
        shared || push!(distinct, x)
    end
    return d
end

function provider(
    reg_path :: AbstractString = REG_PATH;
    excludes :: AbstractSet{<:AbstractString} = EXCLUDES,
)
    reg_inst = RegistryInstance(reg_path)
    reg_dict = Dict(p.name => p
        for p in values(reg_inst.pkgs) if p.name ∉ excludes)

    DepsProvider(keys(reg_dict)) do pkg :: AbstractString
        # Julia 1.14 dropped the one-arg `init_package_info!(::PkgEntry)` in
        # favor of `init_package_info!(::RegistryInstance, ::PkgEntry)`.
        entry = reg_dict[pkg]
        info = applicable(init_package_info!, reg_inst, entry) ?
            init_package_info!(reg_inst, entry) : init_package_info!(entry)
        vers = sort!(collect(keys(info.version_info)), rev=true)
        deps = Dict(v => String[] for v in vers)
        comp = Dict(v => Dict{String,VersionSpec}() for v in vers)
        # scan versions and populate deps & compat data
        for v in vers
            for (r, d) in info.deps
                v in r && union!(deps[v], _dep_names(d, reg_inst))
            end
            for (r, c) in info.compat
                v in r && mergewith!(∩, comp[v], _compat_names(c, reg_inst))
            end
            hasfield(typeof(info), :weak_compat) &&
            for (r, c) in info.weak_compat
                v in r && mergewith!(∩, comp[v], _compat_names(c, reg_inst))
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
        dedup_values!(deps, vers)
        dedup_values!(comp, vers)
        # return resolver PkgData structure
        PkgData(vers, deps, comp)
    end
end

end # module
