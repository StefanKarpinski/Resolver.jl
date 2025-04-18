module Registries

export registry_provider

import Base: UUID
import HistoricalStdlibVersions # populates data for get_last_stdlibs
import JSON
import Pkg.Registry: JULIA_UUID, PkgEntry, init_package_info!, reachable_registries
import Pkg.Types: get_last_stdlibs
import Pkg.Versions: VersionSpec
import Resolver: DepsProvider, PkgData

## load Julia versions

julia_versions_data = JSON.parsefile(julia_versions_file)
const JULIA_VERSIONS = sort!(VersionNumber.(keys(julia_versions_data)))

# drop prerelease versions
# unless there is no corresponding final release
# in which case only keep the last prerelease
filter!(JULIA_VERSIONS) do v
    isempty(v.prerelease) && return true
    Base.thispatch(v) ∈ JULIA_VERSIONS && return false
    all(JULIA_VERSIONS) do v′
        Base.thispatch(v) ≠ Base.thispatch(v′) || v′ ≤ v
    end
end

## extracting the dependency graph from registries

function sort_versions_default(uuid::UUID, vers::Set{VersionNumber})
    sort!(collect(vers), rev=true)
end

function registry_provider(;
    julia_version  :: VersionNumber = VERSION,
    project_compat :: Dict{UUID,VersionSpec} = Dict{UUID,VersionSpec}(),
    sort_versions  :: Function = sort_versions_default,
)
    packages = Dict{UUID,Vector{PkgEntry}}()
    excludes = Set(keys(get_last_stdlibs(julia_version)))

    for reg in reachable_registries()
        for (uuid, entry) in reg.pkgs
            uuid in excludes && continue
            push!(get!(()->PkgEntry[], packages, uuid), entry)
        end
    end

    rp = DepsProvider(keys(packages)) do uuid::UUID
        vers = Set{VersionNumber}()
        deps = Dict{VersionNumber,Vector{UUID}}()
        comp = Dict{VersionNumber,Dict{UUID,VersionSpec}}()
        uuid == JULIA_UUID &&
            return PkgData([julia_version], deps, comp)
        for entry in packages[uuid]
            info = init_package_info!(entry)
            # compat-filtered versions from this registry
            new_vers = if !haskey(project_compat, uuid)
                collect(keys(info.version_info))
            else
                filter(keys(info.version_info)) do v
                    v in project_compat[uuid]
                end
            end
            # scan versions and populate deps & compat data
            for v in new_vers
                push!(vers, v)
                uuids = Dict{String,UUID}()
                deps_v = get!(()->valtype(deps)(), deps, v)
                for (r, d) in info.deps
                    v in r || continue
                    merge!(uuids, d)
                    union!(deps_v, values(d))
                end
                comp_v = get!(()->valtype(comp)(), comp, v)
                for (r, c) in info.compat
                    v in r || continue
                    for (name, spec) in c
                        u = uuids[name]
                        if haskey(comp_v, u)
                            comp_v[u] = spec ∩ comp_v[u]
                        else
                            comp_v[u] = spec
                        end
                    end
                end
            end
        end
        # scrub excluded uuids from deps and sort
        for d in values(deps)
            setdiff!(d, excludes)
            sort!(d)
        end
        # scrub excluded uuids from compat
        for c in values(comp), x in excludes
            delete!(c, x)
        end
        # sort versions
        vers = sort_versions(uuid, vers) :: Vector{VersionNumber}
        # deduplicate data structures to save some memory
        for i = 1:length(vers)-1, j = i+1:length(vers)
            v, w = vers[i], vers[j]
            deps[v] == deps[w] && (deps[v] = deps[w])
            comp[v] == comp[w] && (comp[v] = comp[w])
        end
        # return resolver PkgData structure
        PkgData(vers, deps, comp)
    end

    return packages, rp
end

end # module

using .Registries
