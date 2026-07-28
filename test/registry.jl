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
        # a version's data is determined by which registry entries apply
        # to it, so group versions by their applicable-entry sets and
        # materialize each group's deps & compat once, shared by all of
        # the group's versions (which also deduplicates the structures)
        dep_entries  = collect(info.deps)
        comp_entries = collect(info.compat)
        weak_entries = hasfield(typeof(info), :weak_compat) ?
            collect(info.weak_compat) : empty(comp_entries)
        nd, nc = length(dep_entries), length(comp_entries)
        deps = Dict{VersionNumber,Vector{String}}()
        comp = Dict{VersionNumber,Dict{String,VersionSpec}}()
        groups = Dict{Vector{Int},
                      Tuple{Vector{String},Dict{String,VersionSpec}}}()
        key = Int[]
        for v in vers
            empty!(key)
            for (i, (r, _)) in enumerate(dep_entries)
                v in r && push!(key, i)
            end
            for (i, (r, _)) in enumerate(comp_entries)
                v in r && push!(key, nd + i)
            end
            for (i, (r, _)) in enumerate(weak_entries)
                v in r && push!(key, nd + nc + i)
            end
            group = get(groups, key, nothing)
            if group === nothing
                d = String[]
                c = Dict{String,VersionSpec}()
                for i in key
                    if i ≤ nd
                        union!(d, _dep_names(dep_entries[i][2], reg_inst))
                    elseif i ≤ nd + nc
                        mergewith!(∩, c,
                            _compat_names(comp_entries[i - nd][2], reg_inst))
                    else
                        mergewith!(∩, c,
                            _compat_names(weak_entries[i - nd - nc][2], reg_inst))
                    end
                end
                # sort deps & scrub out excludes (stdlibs, julia itself)
                sort!(d)
                setdiff!(d, excludes)
                for x in excludes
                    delete!(c, x)
                end
                group = (d, c)
                groups[copy(key)] = group
            end
            deps[v], comp[v] = group
        end
        # return resolver PkgData structure
        PkgData(vers, deps, comp)
    end
end

end # module
