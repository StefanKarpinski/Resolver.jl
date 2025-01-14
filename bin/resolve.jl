#!/usr/bin/env julia

## parse command arguments

const USAGE = """
usage: $PROGRAM_FILE [options] [<project path>]

  --julia=<version>       version to resolve for (default: $VERSION)
  --additional=<pkgs>     additional packages to require
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
            additional |
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

## some global constants

const JULIA_VER = !haskey(OPTS, "julia") ? VERSION : let
    # TODO: treat argument as VersionSpec and download versions
    # from https://julialang-s3.julialang.org/bin/versions.json
    # to get the actual Julia versions and pick maximal matching
    str = OPTS["julia"]
    isnothing(str) && usage("Option requires an argument: --julia")
    ver = tryparse(VersionNumber, str)
    @something ver usage("Invalid version number: --julia=$str")
end
const EXCLUDES = Set(keys(get_last_stdlibs(JULIA_VER)))
const PACKAGES = Dict{UUID,Vector{PkgEntry}}()
const UUIDS = Dict{String,Vector{UUID}}()
const COMPAT = Dict{UUID,VersionSpec}() # populated later

for reg in reachable_registries()
    for (uuid, entry) in reg.pkgs
        uuid in EXCLUDES && continue
        push!(get!(()->PkgEntry[], PACKAGES, uuid), entry)
    end
end
for (uuid, entries) in PACKAGES, entry in entries
    push!(get!(()->UUID[], UUIDS, entry.name), uuid)
end
foreach(sort!, values(UUIDS))

## declare functions sorting functions

function sort_versions end
function sort_packages_by end

## extracting the dependency graph from registries

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

let env = EnvCache(PROJ)
    proj = env.project
    global const uuids = merge(proj.deps, proj.weakdeps, proj.extras)
    uuids["julia"] = JULIA_UUID
    global const deps = filter!(∉(EXCLUDES), collect(values(proj.deps)))
    push!(deps, JULIA_UUID)
    for (name, comp) in proj.compat
        COMPAT[uuids[name]] = comp.val
    end
end

## parsing a packages spec

function parse_packages(str::AbstractString)
    pkgs = UUID[]
    for item in split(str, ',')
        if item == "@deps"
            union!(pkgs, deps)
            continue
        end
        uuid = tryparse(UUID, item)
        if isnothing(uuid)
            if item in keys(uuids)
                uuid = uuids[item]
            elseif item in keys(UUIDS)
                list = UUIDS[item]
                length(list) == 1 || usage("Ambiguous package name: $item")
                uuid = only(list)
            end
        end
        uuid isa UUID || usage("Invalid value for --prioritize: $item")
        uuid ∉ pkgs && push!(pkgs, uuid)
    end
    return pkgs
end

## additional requirements

if haskey(OPTS, "additional")
    str = OPTS["additional"]
    isnothing(str) && usage("Option requires an argument: --additional")
    union!(deps, parse_packages(str))
end

## sorting packages (resolution priority)

const ZERO_UUID = reinterpret(UUID, UInt128(0))

function default_sort_packages_by(uuid::UUID)
    uuid == JULIA_UUID ? ZERO_UUID : uuid
end

sort!(deps, by = default_sort_packages_by)

if !haskey(OPTS, "prioritize")
    sort_packages_by(uuid::UUID) = default_sort_packages_by(uuid)
else
    let
        str = OPTS["prioritize"]
        isnothing(str) && usage("Option requires an argument: --prioritize")
        prioritize = parse_packages(str)
        priority = Dict(map(reverse, enumerate(prioritize)))
        ∞ = length(prioritize) + 1
        global function sort_packages_by(uuid::UUID)
            get(priority, uuid, ∞), default_sort_packages_by(uuid)
        end
        @show prioritize
    end
end

## interpret options

for (key, val) in OPTS

end

## sorting versions and packages

function sort_versions(uuid::UUID, vers::Set{VersionNumber})
    sort!(collect(vers), rev=true)
end

## do an actual resolve

pkgs, vers = resolve(dp(), deps; by = sort_packages_by)
names = [first(PACKAGES[u]).name for u in pkgs]
display([names vers])
