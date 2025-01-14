#!/usr/bin/env julia

const USAGE = """
usage: $PROGRAM_FILE [options] [<project path>]

  --julia=<version>       Julia version to resolve for (default: $VERSION)
  --prioritize=<pkgs>     package names/uuids to prioritize
  --fix[=<pkgs>]          prefer current full version number
  --fix-minor[=<pkgs>]    prefer current major.minor version
  --fix-major[=<pkgs>]    prefer current major version
  --max[=<pkgs>]          maximize major.minor.patch
  --max-minor[=<pkgs>]    maximize major.minor; minimize patch
  --max-major[=<pkgs>]    maximize major; minimize minor.patch
  --min[=<pkgs>]          minimize major.minor.patch
  --min-minor[=<pkgs>]    minimize major.minor; maximize patch
  --min-major[=<pkgs>]    minimize major; maximize minor.patch

Wherever <pkgs> appears you can specify a comma separated list of:

  * Package uuids for any packages
  * Package names from any the [deps], [weakdeps] or [extras] sections
    of the specified environment's project file
  * @deps for all packages in the [deps] section of the project file
"""

function usage()
    println(stderr, USAGE)
    exit(0)
end

function usage(msg)
    println(stderr, "[ERROR] $msg\n\n$USAGE")
    exit(1)
end

expand_project(path::Nothing) = Base.active_project()

function expand_project(path::AbstractString)
    isdir(path) && for name in Base.project_names
        file = joinpath(path, name)
        Base.isfile_casesensitive(file) && return file
    end
    Base.isfile_casesensitive(path) && return path
    usage("Not a project path: $path")
end

const OPTS = Dict{String,Union{String,Nothing}}()
const PROJ = let proj = nothing,
    opt_re = r"""
        ^--(
            julia |
            prioritize |
            (?:fix|max|min)(?:-(?:minor|major))?
        )(?:=(.+))?$
    """x
    for arg in ARGS
        if startswith(arg, "-")
            arg in ("-h", "--help") && usage()
            m = match(opt_re, arg)
            isnothing(m) && usage("Invalid option: $arg")
            OPTS[m[1]] = m[2]
        else
            isnothing(proj) || usage("At most one project can be specified.")
            proj = arg
        end
    end
    expand_project(proj)
end

## imports

push!(empty!(LOAD_PATH), @__DIR__)

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

function sort_versions(uuid::UUID, vers::Set{VersionNumber})
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
    vers = sort_versions(uuid, vers)
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

env = EnvCache(PROJ)
proj = env.project
uuids = merge(proj.deps, proj.weakdeps, proj.extras)
uuids["julia"] = JULIA_UUID
deps = sort!(filter!(∉(EXCLUDES), collect(values(proj.deps))))
push!(deps, JULIA_UUID)

# populate COMPAT dict for dependency provider
for (name, comp) in proj.compat
    COMPAT[uuids[name]] = comp.val
end

## do an actual resolve

pkgs, vers = resolve(dp(), deps; by = sort_packages_by)
names = [first(PACKAGES[u]).name for u in pkgs]
display([names vers])
