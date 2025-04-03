#!/usr/bin/env julia

push!(empty!(LOAD_PATH), @__DIR__)

## parsing command arguments

include("Options.jl")

USAGE[] = """
usage: $PROGRAM_FILE [options] [<project path>]

  --help -h               print this help message
  --manifest[=<file>]     write new manifest to file
  --julia=<version>       version to resolve for (default: $VERSION)
  --additional=<pkgs>     additional packages to require
  --prioritize=<pkgs>     package names/uuids to prioritize
  --fix[=<pkgs>]          prefer current full version number
  --fix-minor[=<pkgs>]    prefer current major.minor version
  --fix-major[=<pkgs>]    prefer current major version
  --unfix=[<pkgs>]        undo (override) previous fix options
  --max[=<pkgs>]          maximize version number
  --max-minor[=<pkgs>]    maximize major.minor (minimize patch)
  --max-major[=<pkgs>]    maximize major (minimize minor.patch)
  --min[=<pkgs>]          minimize version number
  --min-minor[=<pkgs>]    minimize major.minor (maximize patch)
  --min-major[=<pkgs>]    minimize major (maximize minor.patch)

Wherever <pkgs> appears you can specify a comma separated list of:

  * Package uuids for any packages
  * Package names from any the [deps], [weakdeps] or [extras] sections
    of the specified environment's project file
  * @deps for all packages in the [deps] section of the project file
"""

parse_opts!(ARGS, split("""
    manifest julia additional prioritize
    fix fix-minor fix-major unfix
    max max-minor max-major
    min min-minor min-major
"""))

length(ARGS) ≤ 1 || usage("At most one project can be specified.")

function expand_project(path::AbstractString)
    isdir(path) && for name in Base.project_names
        file = joinpath(path, name)
        Base.isfile_casesensitive(file) && return file
    end
    Base.isfile_casesensitive(path) && return path
    usage("Not a project path: $path")
end

const PROJ = length(ARGS) ≥ 1 ?
    expand_project(ARGS[1]) :
    Base.active_project()

## imports

include("Registries.jl")

import Base: SHA1, UUID, thismajor, thisminor
import Pkg.Registry:
    JULIA_UUID, PkgEntry, RegistryInstance,
    init_package_info!, reachable_registries
import Pkg.Types:
    Context, EnvCache, Manifest, PackageEntry,
    get_last_stdlibs, write_manifest
import Pkg.Versions: VersionSpec
import Resolver: DepsProvider, PkgData, resolve

## options: target Julia version

const julia_version = handle_opts(:julia, VERSION) do val::String
    # TODO: treat argument as VersionSpec and download versions
    # from https://julialang-s3.julialang.org/bin/versions.json
    # to get the actual Julia versions and pick maximal matching
    ver = tryparse(VersionNumber, val)
    @something ver usage("Invalid version number: --julia=$val")
end

## load project & manifest

const env = EnvCache(PROJ)
const ctx = Context(; env)

let proj = env.project
    global const project_uuids = merge(proj.deps, proj.weakdeps, proj.extras)
    project_uuids["julia"] = JULIA_UUID
    global const project_deps = collect(values(proj.deps))
    push!(project_deps, JULIA_UUID)
    global const project_compat = Dict{UUID,VersionSpec}()
    for (name, comp) in proj.compat
        project_compat[project_uuids[name]] = comp.val
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
            union!(pkgs, project_deps)
            continue
        end
        uuid = tryparse(UUID, item)
        if isnothing(uuid)
            if item in keys(project_uuids)
                uuid = project_uuids[item]
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

global const reqs = copy(project_deps)

handle_opts(:additional) do val::String
    union!(reqs, parse_packages(val))
end

## options: sorting packages (resolution priority)

const ZERO_UUID = UUID(0)

function default_sort_packages_by(uuid::UUID)
    uuid == JULIA_UUID ? ZERO_UUID : uuid
end

sort!(reqs, by = default_sort_packages_by)

sort_packages_by(uuid::UUID) = default_sort_packages_by(uuid)

let i = 0
    priority = Dict{UUID,Int}()
    handle_opts(:prioritize) do str::String
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

# zero UUID used for defaults
const FIX_PATCH = Dict(ZERO_UUID => false)
const FIX_MINOR = Dict(ZERO_UUID => false)
const FIX_MAJOR = Dict(ZERO_UUID => false)
const ORDER_MAP = Dict(ZERO_UUID => :max)

handle_opts(r"^((un)?fix|max|min)") do opt, val
    pkgs = isnothing(val) ? [ZERO_UUID] : parse_packages(val)
    if opt == :unfix
        for uuid in pkgs,
            dict in (FIX_PATCH, FIX_MINOR, FIX_MAJOR)
            dict[uuid] = false
        end
    elseif opt == :fix
        for uuid in pkgs
            FIX_PATCH[uuid] = true
        end
    elseif opt == :fix_minor
        for uuid in pkgs
            FIX_MINOR[uuid] = true
        end
    elseif opt == :fix_major
        for uuid in pkgs
            FIX_MAJOR[uuid] = true
        end
    elseif opt in (:max, :max_minor, :max_major, :min, :min_minor, :min_major)
        for uuid in pkgs
            ORDER_MAP[uuid] = opt
        end
    else
        error("Internal error: unexpected option name: $opt")
    end
end

## version sorting

function fixed(
    u::Nothing, # no old version
    fix_patch::Bool,
    fix_minor::Bool,
    fix_major::Bool,
)
    v::VersionNumber -> 0x0
end

function fixed(
    u::VersionNumber, # old version
    fix_patch::Bool,
    fix_minor::Bool,
    fix_major::Bool,
) :: Function
    function fixed_by(v::VersionNumber)
        fix_patch && u == v                       ? 0x3 :
        fix_minor && thisminor(u) == thisminor(v) ? 0x2 :
        fix_major && thismajor(u) == thismajor(v) ? 0x1 : 0x0
    end
end

function order(level::Symbol) :: Function
    level == :min && return (u::VersionNumber, v::VersionNumber) -> u < v
    level == :max && return (u::VersionNumber, v::VersionNumber) -> u > v
    level == :min_minor && return (u::VersionNumber, v::VersionNumber) ->
        thisminor(u) ≠ thisminor(v) ? thisminor(u) < thisminor(v) : u > v
    level == :max_minor && return (u::VersionNumber, v::VersionNumber) ->
        thisminor(u) ≠ thisminor(v) ? thisminor(u) > thisminor(v) : u < v
    level == :min_major && return (u::VersionNumber, v::VersionNumber) ->
        thismajor(u) ≠ thismajor(v) ? thismajor(u) < thismajor(v) : u > v
    level == :max_major && return (u::VersionNumber, v::VersionNumber) ->
        thismajor(u) ≠ thismajor(v) ? thismajor(u) > thismajor(v) : u < v
end

function sort_versions(uuid::UUID, vers::Set{VersionNumber})
    fixed_by = fixed(
        get(old_versions, uuid, nothing),
        get(FIX_PATCH, uuid, FIX_PATCH[ZERO_UUID]),
        get(FIX_MINOR, uuid, FIX_MINOR[ZERO_UUID]),
        get(FIX_MAJOR, uuid, FIX_MAJOR[ZERO_UUID]),
    )
    order_lt = order(get(ORDER_MAP, uuid, ORDER_MAP[ZERO_UUID]))
    function lt(u::VersionNumber, v::VersionNumber)
        fixed_u = fixed_by(u) :: UInt8
        fixed_v = fixed_by(v) :: UInt8
        fixed_u ≠ fixed_v ? fixed_u > fixed_v : order_lt(u, v)
    end
    sort!(collect(vers); lt)
end

## do an actual resolve

const packages, rp =
    registry_provider(; julia_version, project_compat, sort_versions)
intersect!(reqs, rp.packages)
pkgs, vers = resolve(rp, reqs; max=1, by=sort_packages_by)
for uuid in reqs
    i = findfirst(==(uuid), pkgs)
    isnothing(vers[i]) && error("Unsatisfiable")
end

handle_opts(:manifest, false) do val
    manifest_file = something(val, env.manifest_file)
    # generate a manifest file
    manifest_deps = Dict{UUID,PackageEntry}()
    for (i, uuid) in enumerate(pkgs)
        uuid === JULIA_UUID && continue
        version = vers[i]
        version === nothing && continue
        infos = Set{Tuple{String,SHA1,Dict{String,UUID},Dict{String,UUID}}}()
        for entry in packages[uuid]
            info = init_package_info!(entry)
            haskey(info.version_info, version) || continue
            tree = info.version_info[version].git_tree_sha1
            deps = Dict{String,UUID}()
            for (r, d) in info.deps
                version in r || continue
                merge!(deps, d)
            end
            if get(deps, "julia", nothing) == JULIA_UUID
                delete!(deps, "julia")
            end
            weakdeps = Dict{String,UUID}()
            for (r, d) in info.weak_deps
                version in r || continue
                merge!(weakdeps, d)
            end
            push!(infos, (entry.name, tree, deps, weakdeps))
        end
        if length(infos) ≠ 1
            name = first(packages[uuid]).name
            abbr = string(uuid)[1:8]
            error("Package $name [$abbr]: version $version resolved but " *
                (length(infos) > 1 ? "has conflicting definitions" : "not found"))
        end
        name, tree_hash, deps, weakdeps = only(infos)
        manifest_deps[uuid] = PackageEntry(;
            name,
            version,
            tree_hash,
            deps,
            weakdeps,
        )
    end
    # TODO: need to handle project_hash, stdlibs, deps, weakdeps & extensions
    manifest = Manifest(; julia_version, deps = manifest_deps)
    if manifest_file == "-"
        write_manifest(stdout, manifest)
    else
        write_manifest(manifest, manifest_file)
    end
    return true
end || begin
    # just print packages and versions
    names = [first(packages[uuid]).name for uuid in pkgs]
    width = maximum(length, names)
    for (i, uuid) in enumerate(pkgs)
        uuid === JULIA_UUID && continue
        version = vers[i]
        version === nothing && continue
        name = first(packages[uuid]).name
        try println(uuid, " ", rpad(name, width), " ", version)
        catch err
            # no stack trace for SIGPIPE
            err isa Base.IOError && err.code == Base.UV_EPIPE && exit(2)
            rethrow() # some other error
        end
    end
end
