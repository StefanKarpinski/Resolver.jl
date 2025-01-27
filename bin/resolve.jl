#!/usr/bin/env julia

## parsing command arguments

const USAGE = """
usage: $PROGRAM_FILE [options] [<project path>]

  --save[=<file>]         write new manifest to file (default: in-place)
  --julia=<version>       version to resolve for (default: $VERSION)
  --additional=<pkgs>     additional packages to require
  --prioritize=<pkgs>     package names/uuids to prioritize
  --fix[=<pkgs>]          prefer current full version number
  --fix-minor[=<pkgs>]    prefer current major.minor version
  --fix-major[=<pkgs>]    prefer current major version
  --unfix=[<pkgs>]        override previous fix options
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

const OPTS = Vector{Pair{String,Union{String,Nothing}}}()
const PROJ = let proj = nothing,
    opt_re = r"""^--(
        save | julia | additional | prioritize |
        unfix | (?:fix|max|min) (?:-(?:minor|major))?
    )(?:=(.+))?$"""x
    for arg in ARGS
        if startswith(arg, "-")
            arg in ("-h", "--help") && usage()
            m = match(opt_re, arg)
            isnothing(m) && usage("Invalid option: $arg")
            push!(OPTS, m[1] => m[2])
        else
            isnothing(proj) || usage("At most one project can be specified.")
            proj = arg
        end
    end
    expand_project(proj)
end

function handle_opt(body::Function, name::AbstractString, value::Any=nothing)
    for (key, val) in OPTS
        key == name || continue
        value = body(val)
    end
    return value
end

## imports

push!(empty!(LOAD_PATH), @__DIR__)

import Base: SHA1, UUID, thismajor, thisminor
import HistoricalStdlibVersions
import Pkg.Registry:
    JULIA_UUID, PkgEntry, RegistryInstance,
    init_package_info!, reachable_registries
import Pkg.Types:
    EnvCache, Manifest, PackageEntry,
    get_last_stdlibs, write_manifest
import Pkg.Versions: VersionSpec
import Resolver: DepsProvider, PkgData, resolve

## options: target Julia version

const julia_version = handle_opt("julia", VERSION) do str
    isnothing(str) && usage("Option requires an argument: --julia")
    # TODO: treat argument as VersionSpec and download versions
    # from https://julialang-s3.julialang.org/bin/versions.json
    # to get the actual Julia versions and pick maximal matching
    ver = tryparse(VersionNumber, str)
    @something ver usage("Invalid version number: --julia=$str")
end

const EXCLUDES = Set(keys(get_last_stdlibs(julia_version)))
const PACKAGES = Dict{UUID,Vector{PkgEntry}}()
const UUIDS = Dict{String,Vector{UUID}}()

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

## load project & manifest

const env = EnvCache(PROJ)
let proj = env.project
    global const uuids = merge(proj.deps, proj.weakdeps, proj.extras)
    uuids["julia"] = JULIA_UUID
    global const deps = filter!(∉(EXCLUDES), collect(values(proj.deps)))
    push!(deps, JULIA_UUID)
    global const COMPAT = Dict{UUID,VersionSpec}()
    for (name, comp) in proj.compat
        COMPAT[uuids[name]] = comp.val
    end
    global const old_versions = Dict{UUID,VersionNumber}()
    for (uuid, entry) in env.manifest.deps
        isnothing(entry.version) && continue
        old_versions[uuid] = entry.version
    end
end

## options: parsing packages specs

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

## options: additional requirements

handle_opt("additional") do str
    isnothing(str) && usage("Option requires an argument: --additional")
    union!(deps, parse_packages(str))
end

## options: sorting packages (resolution priority)

const ZERO_UUID = reinterpret(UUID, UInt128(0))

function default_sort_packages_by(uuid::UUID)
    uuid == JULIA_UUID ? ZERO_UUID : uuid
end

sort!(deps, by = default_sort_packages_by)

sort_packages_by(uuid::UUID) = default_sort_packages_by(uuid)

let i = 0
    priority = Dict{UUID,Int}()
    handle_opt("prioritize") do str
        isnothing(str) && usage("Option requires an argument: --prioritize")
        for uuid in parse_packages(str)
            priority[uuid] = (i += 1)
        end
    end
    if !isempty(priority)
        n = i + 1
        global function sort_packages_by(uuid::UUID)
            get(priority, uuid, n), default_sort_packages_by(uuid)
        end
    end
end

## options: sorting versions

FIXED_DEF::String = "none"
ORDER_DEF::String = "max"

const FIXED_MAP = Dict{UUID,String}()
const ORDER_MAP = Dict{UUID,String}()

fix_level(key::AbstractString) =
    key == "unfix" ? "none" :
    key == "fix"   ? "all" :
    chopprefix(key, "fix-")

for (key, val) in OPTS
    startswith(key, r"(un)?fix|max|min") || continue
    if isnothing(val)
        # set defaults
        if contains(key, "fix")
            global FIXED_DEF = fix_level(key)
        else
            global ORDER_DEF = key
        end
    else
        pkgs = parse_packages(val)
        for uuid in pkgs
            if contains(key, "fix")
                FIXED_MAP[uuid] = fix_level(key)
            else
                ORDER_MAP[uuid] = key
            end
        end
    end
end

## helpers for version sorting

fixed(level::String, u::Nothing) = v::VersionNumber -> true

function fixed(level::String, u::VersionNumber) :: Function
    level == "none"  && return v::VersionNumber -> true
    level == "all"   && return v::VersionNumber -> v == u
    level == "minor" && return v::VersionNumber -> thisminor(v) == thisminor(u)
    level == "major" && return v::VersionNumber -> thismajor(v) == thismajor(u)
end

function order(level::String) :: Function
    level == "min" && return (u::VersionNumber, v::VersionNumber) -> u < v
    level == "max" && return (u::VersionNumber, v::VersionNumber) -> u > v
    level == "min-minor" && return (u::VersionNumber, v::VersionNumber) ->
        thisminor(u) ≠ thisminor(v) ? thisminor(u) < thisminor(v) : u > v
    level == "max-minor" && return (u::VersionNumber, v::VersionNumber) ->
        thisminor(u) ≠ thisminor(v) ? thisminor(u) > thisminor(v) : u < v
    level == "min-major" && return (u::VersionNumber, v::VersionNumber) ->
        thismajor(u) ≠ thismajor(v) ? thismajor(u) < thismajor(v) : u > v
    level == "max-major" && return (u::VersionNumber, v::VersionNumber) ->
        thismajor(u) ≠ thismajor(v) ? thismajor(u) > thismajor(v) : u < v
end

## extracting the dependency graph from registries

dp() = DepsProvider(keys(PACKAGES)) do uuid::UUID
    vers = Set{VersionNumber}()
    deps = Dict{VersionNumber,Vector{UUID}}()
    comp = Dict{VersionNumber,Dict{UUID,VersionSpec}}()
    uuid == JULIA_UUID &&
        return PkgData([julia_version], deps, comp)
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
    fixed_p = fixed(
        get(FIXED_MAP, uuid, FIXED_DEF),
        get(old_versions, uuid, nothing),
    )
    order_lt = order(get(ORDER_MAP, uuid, ORDER_DEF))
    function lt(u::VersionNumber, v::VersionNumber)
        fixed_u = fixed_p(u)
        fixed_v = fixed_p(v)
        fixed_u ≠ fixed_v ? fixed_u > fixed_v : order_lt(u, v)
    end
    vers = sort!(collect(vers); lt)
    # deduplicate data structures to save some memory
    for i = 1:length(vers)-1, j = i+1:length(vers)
        v, w = vers[i], vers[j]
        deps[v] == deps[w] && (deps[v] = deps[w])
        comp[v] == comp[w] && (comp[v] = comp[w])
    end
    # return resolver PkgData structure
    PkgData(vers, deps, comp)
end

## do an actual resolve

versions = let
    pkgs, vers = resolve(dp(), deps; max=1, by=sort_packages_by)
    for uuid in deps
        i = findfirst(==(uuid), pkgs)
        isnothing(vers[i]) && error("Unsatisfiable")
    end
    Dict{UUID,VersionNumber}(p => v for (p, v) in zip(pkgs, vers))
end

## generate a manifest to output

manifest_deps = Dict{UUID, PackageEntry}()
for (uuid, version) in versions
    uuid == JULIA_UUID && continue
    infos = Set{Tuple{String,SHA1}}()
    for entry in PACKAGES[uuid]
        info = init_package_info!(entry)
        haskey(info.version_info, version) || continue
        tree = info.version_info[version].git_tree_sha1
        push!(infos, (entry.name, tree))
    end
    if length(infos) ≠ 1
        name = first(PACKAGES[uuid]).name
        abbr = string(uuid)[1:8]
        error("Package $name [$abbr]: version $version resolved but " *
            length(infos) > 1 ? "has conflicting definitions" : "not found")
    end
    name, tree_hash = only(infos)
    manifest_deps[uuid] = PackageEntry(; name, version, tree_hash)
end

# TODO: need to handle project_hash, stdlibs and deps graph
manifest = Manifest(; julia_version, deps = manifest_deps)

handle_opt("save", false) do str
    manifest_file = something(str, env.manifest_file)
    write_manifest(manifest, manifest_file)
    return true
end || begin
    write_manifest(stdout, manifest)
end
