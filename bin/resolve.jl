#!/usr/bin/env julia +1.11

push!(empty!(LOAD_PATH), @__DIR__)

## parsing command arguments

include("Options.jl")

USAGE[] = """
usage: $PROGRAM_FILE [options] [<project path>]

  --help -h               print this help message

  --print-manifest        print the new manifest to stdout
  --print-versions        print the resolved versions to stdout

  --julia=<versions>      Julia versions to resolve for (default: 1+)
                          use registry compat syntax, not semver

  --allow-pre[=<pkgs>]    allow prerelease versions
  --extra-deps=<pkgs>     extra packages to require
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
    print-manifest print-versions
    julia allow-pre extra-deps prioritize
    fix fix-minor fix-major unfix
    max max-minor max-major
    min min-minor min-major
"""))

length(ARGS) ≤ 1 || usage("At most one project can be specified.")

output = handle_opts(
    opt::Symbol -> opt,
    r"^print_(manifest|versions)$",
    :write_manifest,
)

function expand_project(path::AbstractString)
    isdir(path) && for name in Base.project_names
        file = joinpath(path, name)
        Base.isfile_casesensitive(file) && return file
    end
    Base.isfile_casesensitive(path) && return path
    usage("Not a project path: $path")
end

const PROJ = length(ARGS) ≥ 1 ? expand_project(ARGS[1]) : Base.active_project()

## imports

include("Registries.jl")

import Base: SHA1, UUID, thismajor, thisminor
import Pkg
import Pkg.Operations: record_project_hash, download_source
if isdefined(Pkg.Operations, :fixups_from_projectfile!)
    import Pkg.Operations: fixups_from_projectfile!
elseif isdefined(Pkg.Operations, :fixup_ext!)
    import Pkg.Operations: fixup_ext!
else
    error("Pkg too old to support generating manifests with extensions")
end
import Pkg.Registry: JULIA_UUID, PkgEntry, RegistryInstance,    init_package_info!, reachable_registries
import Pkg.Types: Context, EnvCache, Manifest, PackageEntry, get_last_stdlibs, write_manifest
import Pkg.Versions: VersionSpec, semver_spec
import Resolver: DepsProvider, PkgData, resolve
import HistoricalStdlibVersions: STDLIBS_BY_VERSION, UNREGISTERED_STDLIBS
import TOML

## options: target Julia version

const julia_versions = handle_opts(:julia, VersionSpec("1")) do val::String
    try VersionSpec(split(val, r"\s*,\s*"))
    catch
        usage("Invalid compat version spec: --julia=$val")
    end
end

## load project & manifest

const env = EnvCache(PROJ)

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
    if env.manifest.julia_version !== nothing
        old_versions[JULIA_UUID] = env.manifest.julia_version
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

handle_opts(:extra_deps) do val::String
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
const allow_pre = Dict(ZERO_UUID => false)
const FIX_PATCH = Dict(ZERO_UUID => false)
const FIX_MINOR = Dict(ZERO_UUID => false)
const FIX_MAJOR = Dict(ZERO_UUID => false)
const ORDER_MAP = Dict(ZERO_UUID => :max)

handle_opts(r"^(allow_pre|(un)?fix|max|min)") do opt, val
    pkgs = isnothing(val) ? [ZERO_UUID] : parse_packages(val)
    if opt == :allow_pre
        for uuid in pkgs
            allow_pre[uuid] = true
        end
    elseif opt == :unfix
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
    <ₒ = order(get(ORDER_MAP, uuid, ORDER_MAP[ZERO_UUID]))
    function lt(u::VersionNumber, v::VersionNumber)
        fixed_u = fixed_by(u) :: UInt8
        fixed_v = fixed_by(v) :: UInt8
        fixed_u ≠ fixed_v ? fixed_u > fixed_v : u <ₒ v
    end
    sort!(collect(vers); lt)
end

## do an actual resolve

const packages, rp =
    registry_provider(; julia_versions, project_compat, sort_versions, allow_pre)
intersect!(reqs, rp.packages)
pkgs, vers = resolve(rp, reqs; max=1, by=sort_packages_by)
for uuid in reqs
    i = findfirst(==(uuid), pkgs)
    isnothing(vers[i]) && error("Unsatisfiable")
end

const julia_version = vers[findfirst(==(JULIA_UUID), pkgs)]
const stdlibs = let last_stdlibs = UNREGISTERED_STDLIBS
    for (v, this_stdlibs) in STDLIBS_BY_VERSION
        v > Base.thispatch(julia_version) && break
        last_stdlibs = this_stdlibs
    end
    last_stdlibs
end

struct ManifestEntry
    name :: String
    version :: Union{Nothing,VersionNumber}
    tree_hash :: Union{Nothing,SHA1}
    deps :: Dict{String,UUID}
    weakdeps :: Dict{String,UUID}
end

const info_map = Dict{UUID,ManifestEntry}()

for (i, uuid) in enumerate(pkgs)
    uuid === JULIA_UUID && continue
    if uuid in keys(stdlibs)
        info = stdlibs[uuid]
        deps = Dict{String,UUID}()
        for dep in info.deps
            dep == JULIA_UUID && continue
            name = stdlibs[dep].name
            deps[name] = dep
        end
        weakdeps = Dict{String,UUID}()
        for dep in info.weakdeps
            dep == JULIA_UUID && continue
            name = stdlibs[dep].name
            weakdeps[name] = dep
        end
        info_map[uuid] = ManifestEntry(
            info.name,
            info.version, # can be nothing
            nothing, # tree hash must be nothing
            deps,
            weakdeps,
        )
    elseif uuid in keys(packages)
        version = vers[i]
        infos = Set{ManifestEntry}()
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
            push!(infos, ManifestEntry(
                entry.name,
                version,
                tree,
                deps,
                weakdeps,
            ))
        end
        if length(infos) < 1
            error("Package $uuid: version $version resolved but no registry entries found")
        elseif length(infos) > 1
            names = unique(info.name for info in infos)
            if length(names) == 1
                name = only(names)
            else
                name = join(sort!(names), "/")
            end
            error("Package $name [$uuid]: version $version resolved but has multiple conflicting definitions")
        end
        info_map[uuid] = only(infos)
    else
        error("Package $uuid resolved but not found (shouldn't happen)")
    end
end

if output == :print_versions
    # just print packages and versions
    width = maximum(textwidth(info.name) for info in values(info_map))
    for (i, uuid) in enumerate(pkgs)
        if uuid == JULIA_UUID
            name = "julia"
            version = julia_version
        else
            name = info_map[uuid].name
            version = something(info_map[uuid].version, julia_version)
        end
        try println(uuid, " ", rpad(name, width), " ", version)
        catch err
            # no stack trace for SIGPIPE
            err isa Base.IOError && err.code == Base.UV_EPIPE && exit(2)
            rethrow() # some other error
        end
    end
else # generate a manifest
    deps = Dict{UUID,PackageEntry}(
        uuid => PackageEntry(;
            uuid,
            info.name,
            info.version,
            info.tree_hash,
            info.deps,
            info.weakdeps,
        )
        for (uuid, info) in info_map
    )
    # create manifest and record project hash
    manifest = Manifest(; julia_version, deps)
    env.manifest = manifest
    record_project_hash(env)
    if julia_version < v"1.6.2"
        manifest.manifest_format = v"1"
    end
    if julia_version ≥ v"1.9"
        # getting extension info requires downloading packages
        # this half-installs packages, so don't pollute the real depot
        push!(DEPOT_PATH, mktempdir())
        ctx = Context(; env)
        download_source(ctx)
        # metaprogram around Pkg internal API differences for
        # injecting extension information into the manifest
        if @isdefined fixups_from_projectfile!
            if applicable(fixups_from_projectfile!, ctx)
                fixups_from_projectfile!(ctx)
            elseif applicable(fixups_from_projectfile!, env)
                fixups_from_projectfile!(env)
            else
                error("Pkg: don't know how to call fixups_from_projectfile!")
            end
        elseif @isdefined fixup_ext!
            if applicable(fixup_ext!, env)
                fixup_ext!(env)
            elseif applicable(fixup_ext!, env, values(manifest))
                fixup_ext!(env, values(manifest))
            else
                error("Pkg: don't know how to call fixup_ext!")
            end
        else
            error("Pkg too old to support generating manifests with extensions")
        end
        pop!(DEPOT_PATH)
    end
    # output the manifest
    if manifest.manifest_format ≥ v"2"
        if output == :print_manifest
            write_manifest(stdout, manifest)
        elseif output == :write_manifest
            write_manifest(env)
        else
            error("internal error: unexpected output format: $output")
        end
    else # avoid warning and include comment with Julia version
        if output == :print_manifest
            io = stdout
        elseif output == :write_manifest
            io = open(env.manifest_file, write=true)
        else
            error("internal error: unexpected output format: $output")
        end
        header = """
        # This file is machine-generated - editing it directly is not advised
        #
        # julia_version = "$julia_version"
        # manifest_format = "$(manifest.manifest_format)"
        # project_hash = "$(env.manifest.other["project_hash"])"
        #

        """
        print(io, header)
        raw_manifest = Pkg.Types.destructure(manifest)
        TOML.print(io, raw_manifest, sorted=true) do x
            (typeof(x) in [String, Nothing, UUID, SHA1, VersionNumber]) && return string(x)
            error("unhandled type `$(typeof(x))`")
        end
        io !== stdout && close(io)
    end
end
