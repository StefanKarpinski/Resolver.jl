#!/usr/bin/env julia

push!(empty!(LOAD_PATH), @__DIR__)

## imports

import Base: UUID
import HistoricalStdlibVersions
import Pkg.Registry:
    init_package_info!,
    JULIA_UUID,
    PkgEntry,
    reachable_registries,
    RegistryInstance
import Pkg.Types: EnvCache, get_last_stdlibs
import Pkg.Versions: VersionSpec
import Resolver: DepsProvider, PkgData, resolve

## sorting of versions and packages

function sort_versions(vers::Set{VersionNumber})
    sort!(collect(vers), rev=true)
end

const ZERO_UUID = reinterpret(UUID, UInt128(0))

function sort_packages_by(uuid::UUID)
    uuid == JULIA_UUID ? ZERO_UUID : uuid
end

## extracting the dependency graph from registries

const JULIA_VER = VERSION # Julia version desired
const EXCLUDES = Set(keys(get_last_stdlibs(JULIA_VER)))
const PACKAGES = Dict{UUID,Vector{PkgEntry}}()
const COMPAT = Dict{UUID,VersionSpec}() # populated later

for reg in reachable_registries()
    for (uuid, entry) in reg.pkgs
        uuid in EXCLUDES && continue
        push!(get!(()->PkgEntry[], PACKAGES, uuid), entry)
    end
end

dp() = DepsProvider(keys(PACKAGES)) do uuid::UUID
    vers = Set{VersionNumber}()
    deps = Dict{VersionNumber,Vector{UUID}}()
    comp = Dict{VersionNumber,Dict{UUID,VersionSpec}}()
    uuid == JULIA_UUID &&
        return PkgData([JULIA_VER], deps, comp)
    for entry in PACKAGES[uuid]
        info = init_package_info!(entry)
        # compat-filtered versions from this registry
        new_vers = if !haskey(COMPAT, uuid)
            collect(keys(info.version_info))
        else
            filter(keys(info.version_info)) do v
                v in COMPAT[uuid]
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

## load project & manifest

path = expanduser("~/dev/HTTP/Project.toml")
env = EnvCache(path)
proj = env.project
uuids = merge(proj.deps, proj.weakdeps, proj.extras)
uuids["julia"] = JULIA_UUID
reqs = sort!(filter!(∉(EXCLUDES), collect(values(proj.deps))))
push!(reqs, JULIA_UUID)

# populate COMPAT dict for dependency provider
for (name, comp) in proj.compat
    COMPAT[uuids[name]] = comp.val
end

## do an actual resolve

pkgs, vers = resolve(dp(), reqs; by = sort_packages_by)
names = [first(PACKAGES[u]).name for u in pkgs]
display([names vers])
