#!/usr/bin/env julia

import Base: UUID
import Pkg.Registry:
    init_package_info!,
    JULIA_UUID,
    PkgEntry,
    reachable_registries,
    RegistryInstance
import Pkg.Types: stdlibs
import Pkg.Versions: VersionSpec
import Resolver: DepsProvider, PkgData, resolve

function sort_versions(vers::Set{VersionNumber})
    sort!(collect(vers), rev=true)
end

const EXCLUDES = push!(Set(keys(stdlibs())), JULIA_UUID)
const PACKAGES = Dict{UUID,Vector{PkgEntry}}()

for reg in reachable_registries()
    for (uuid, entry) in reg.pkgs
        uuid in EXCLUDES && continue
        push!(get!(()->PkgEntry[], PACKAGES, uuid), entry)
    end
end

const dp = DepsProvider(keys(PACKAGES)) do uuid::UUID
    vers = Set{VersionNumber}()
    deps = Dict{VersionNumber,Vector{UUID}}()
    comp = Dict{VersionNumber,Dict{UUID,VersionSpec}}()
    for entry in PACKAGES[uuid]
        info = init_package_info!(entry)
        # add versions from this registry
        union!(vers, keys(info.version_info))
        # scan versions and populate deps & compat data
        for v in keys(info.version_info)
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
        setdiff!(d, EXCLUDES)
        sort!(d)
    end
    # scrub excluded uuids from compat
    for c in values(comp), x in EXCLUDES
        delete!(c, x)
    end
    # sort versions
    vers = sort_versions(vers)
    # deduplicate data structures to save some memory
    for i = 1:length(vers)-1, j = i+1:length(vers)
        v, w = vers[i], vers[j]
        deps[v] == deps[w] && (deps[v] = deps[w])
        comp[v] == comp[w] && (comp[v] = comp[w])
    end
    # return resolver PkgData structure
    PkgData(vers, deps, comp)
end

reqs = [UUID("a93c6f00-e57d-5684-b7b6-d8193f3e46c0")] # DataFrames
pkgs, vers = resolve(dp, reqs)
names = [first(PACKAGES[u]).name for u in pkgs]
