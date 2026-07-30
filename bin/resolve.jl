#!/usr/bin/env julia +1.11

push!(empty!(LOAD_PATH), @__DIR__)

## parsing command arguments

include("Options.jl")

USAGE[] = """
usage: $PROGRAM_FILE [options] [<project path>]

  --help -h               print this help message

  --print-manifest        print the new manifest to stdout
  --print-versions        print the resolved versions to stdout

  --julia=<versions>      Julia versions to resolve for; overrides the
                          project's own [compat] julia bound, which is
                          the default (and 1+ if it has none too)
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

Wherever <pkgs> appears you can specify a comma-separated list of:

  * Package uuids for any packages
  * Package names from any the [deps], [weakdeps] or [extras] sections
    of the specified environment's project file
  * Any package name that is associated with a unique UUID across the
    registries that are currently installed
  * @deps for packages in the [deps] section of Project.toml
  * @weakdeps for packages in the [weakdeps] section of Project.toml
  * @alldeps for packages in both the [deps] and [weakdeps] sections
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
import Resolver: Resolver, DepsProvider, PkgData, Problem, resolve
import HistoricalStdlibVersions: STDLIBS_BY_VERSION, UNREGISTERED_STDLIBS
import TOML

## load project & manifest

const env = EnvCache(PROJ)

# Julia 1.12+ resolves all projects in a workspace into one shared manifest;
# older Pkg has no `workspace` field on EnvCache
const in_workspace = hasfield(typeof(env), :workspace) && !isempty(env.workspace)

# uuids of all workspace projects that are packages (the root and any member
# with a uuid); these are provided by the workspace itself as path entries in
# the manifest and must never be resolved from the registries -- not even
# when one of them shadows a registered package (a local dev copy is fixed
# at its local version, exactly as Pkg treats it)
const workspace_uuids = Set{UUID}()
if in_workspace
    env.project.uuid !== nothing && push!(workspace_uuids, env.project.uuid)
    for (_, ws_proj) in env.workspace
        ws_proj.uuid !== nothing && push!(workspace_uuids, ws_proj.uuid)
    end
end

# fixed resolver nodes for the workspace packages: uuid => (name, local
# version, resolvable dependency uuids). Passed to registry_provider so a
# registry package that depends on a workspace package -- most importantly
# one that shadows a registered package -- resolves against the fixed local
# copy; also used to print versions for workspace packages.
const workspace_pkgs = Dict{UUID,Tuple{String,VersionNumber,Vector{UUID}}}()

let proj = env.project
    global const project_names = merge(proj.deps, proj.weakdeps, proj.extras)
    project_names["julia"] = JULIA_UUID
    global const project_deps = collect(values(proj.deps))
    push!(project_deps, JULIA_UUID)
    global const project_weakdeps = collect(values(proj.weakdeps))
    global const project_compat = Dict{UUID,VersionSpec}()
    for (name, comp) in proj.compat
        project_compat[project_names[name]] = comp.val
    end
    global const old_versions = Dict{UUID,VersionNumber}()
    for (uuid, entry) in env.manifest.deps
        isnothing(entry.version) && continue
        old_versions[uuid] = entry.version
    end
    if env.manifest.julia_version !== nothing
        old_versions[JULIA_UUID] = env.manifest.julia_version
    elseif isfile(env.manifest_file)
        for line in eachline(env.manifest_file)
            m = match(r"^#\s+julia_version\s*=\s*\"(.*?)\"\s*$", line)
            m === nothing && continue
            v = tryparse(VersionNumber, m[1])
            v === nothing && continue
            old_versions[JULIA_UUID] = v
        end
    end
end

## load packages from installed registries

const packages = Dict{UUID,Vector{PkgEntry}}()

for reg in reachable_registries()
    for (uuid, entry) in reg.pkgs
        push!(get!(()->PkgEntry[], packages, uuid), entry)
    end
end

const package_names = Dict{String,Vector{UUID}}()

for (uuid, entries) in packages, entry in entries
    uuids = get!(()->valtype(package_names)(), package_names, entry.name)
    uuid in uuids || push!(uuids, uuid)
end
foreach(sort!, values(package_names))

## load workspace sub-project dependencies

if in_workspace
    # collect all resolvable UUIDs: registries + stdlibs, excluding the
    # workspace's own packages
    let known = Set{UUID}(keys(packages))
        push!(known, JULIA_UUID)
        for (_, this_stdlibs) in STDLIBS_BY_VERSION
            union!(known, keys(this_stdlibs))
        end
        union!(known, keys(UNREGISTERED_STDLIBS))
        setdiff!(known, workspace_uuids)
        # record the fixed node for every workspace package: local version
        # plus the local deps that are resolvable or workspace-internal
        let ws_projs = Any[env.project]
            for (_, ws_proj) in env.workspace
                push!(ws_projs, ws_proj)
            end
            for proj in ws_projs
                proj.uuid === nothing && continue
                node_deps = UUID[u for u in values(proj.deps)
                                 if u in known || u in workspace_uuids]
                sort!(node_deps)
                workspace_pkgs[proj.uuid] = (proj.name,
                    something(proj.version, v"0.0.0"), node_deps)
            end
        end
        # the root's deps/weakdeps/compat were collected above without knowing
        # about the workspace; drop any that refer to workspace packages (they
        # are satisfied by the workspace's path entries, not by resolution)
        filter!(∉(workspace_uuids), project_deps)
        filter!(∉(workspace_uuids), project_weakdeps)
        # compat on a workspace package must still admit its fixed local
        # version, as Pkg enforces for fixed packages -- but it must not
        # constrain a same-uuid registry package, so drop it after checking
        for uuid in workspace_uuids
            spec = pop!(project_compat, uuid, nothing)
            spec === nothing && continue
            name, v, _ = workspace_pkgs[uuid]
            v in spec || error(
                "compat with workspace package $name = \"$spec\" excludes its local version $v")
        end
        for (_, ws_proj) in env.workspace
            # add name → UUID mappings for reference (if two workspace projects
            # map the same name to different uuids, the last one wins and any
            # compat on that name is misattributed; Pkg-side collisions like
            # that are pathological, so we don't try harder)
            merge!(project_names, ws_proj.deps, ws_proj.weakdeps, ws_proj.extras)
            # add resolvable deps
            for uuid in values(ws_proj.deps)
                uuid in known && uuid ∉ project_deps && push!(project_deps, uuid)
            end
            # add resolvable weakdeps
            for uuid in values(ws_proj.weakdeps)
                uuid in known && uuid ∉ project_weakdeps && push!(project_weakdeps, uuid)
            end
            # intersect compat constraints; entries naming a workspace package
            # describe the local copy -- check they admit its fixed local
            # version (as Pkg does) but don't let them constrain a same-uuid
            # registry package
            for (name, comp) in ws_proj.compat
                haskey(project_names, name) || continue
                uuid = project_names[name]
                if uuid in workspace_uuids
                    _, v, _ = workspace_pkgs[uuid]
                    v in comp.val || error(
                        "compat with workspace package $name = \"$(comp.val)\" excludes its local version $v")
                    continue
                end
                if haskey(project_compat, uuid)
                    project_compat[uuid] = intersect(project_compat[uuid], comp.val)
                else
                    project_compat[uuid] = comp.val
                end
            end
        end
    end
end

## options: target Julia version
#
# The Julia versions to resolve for come from `--julia` if given, otherwise
# from the project's own `[compat] julia` bound (intersected across the
# workspace above), otherwise from the `1` default. `--julia` overrides the
# project bound outright rather than intersecting with it, so the flag always
# says exactly what will be resolved for.
#
# Whichever source wins feeds the one `julia_versions` spec that determines
# both the `julia` package's version universe and, through it, which
# historical stdlib versions are bundled and pinned (see Registries.jl) --
# so a project bound and the equivalent `--julia` are indistinguishable
# downstream. The bound is *consumed* here: it must not also constrain
# `julia` as a user compat entry in the `Problem` below, or the two readings
# of the same entry would compound.

const julia_versions = something(
    handle_opts(:julia) do val::String
        try VersionSpec(split(val, r"\s*,\s*"))
        catch
            usage("Invalid compat version spec: --julia=$val")
        end
    end,
    pop!(project_compat, JULIA_UUID, nothing),
    VersionSpec("1"),
)

## options: parsing packages specs

function parse_packages(str::AbstractString)
    pkgs = UUID[]
    for item in split(str, ',')
        if item in ("@deps", "@alldeps")
            union!(pkgs, project_deps)
            continue
        end
        if item in ("@weakdeps", "@alldeps")
            union!(pkgs, project_weakdeps)
            continue
        end
        uuid = tryparse(UUID, item)
        if isnothing(uuid) # must be a name then
            if item in keys(project_names)
                uuid = project_names[item]
            elseif item in keys(package_names)
                uuids = package_names[item]
                length(uuids) == 1 || usage("Ambiguous package name: $item")
                uuid = only(uuids)
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

## options: sorting versions

const ZERO_UUID = UUID(0)

# zero UUID used for defaults
const allow_pre = Dict(ZERO_UUID => false)
const FIX_PATCH = Dict(ZERO_UUID => false)
const FIX_MINOR = Dict(ZERO_UUID => false)
const FIX_MAJOR = Dict(ZERO_UUID => false)
const ORDER_MAP = Dict(ZERO_UUID => :max)

# list of all fixed packages
const FIXED = UUID[]

handle_opts(r"^(allow_pre|(un)?fix|max|min)") do opt, val
    pkgs = isnothing(val) ? [ZERO_UUID] : parse_packages(val)
    if opt == :allow_pre
        for uuid in pkgs
            allow_pre[uuid] = true
        end
    elseif opt == :fix
        for uuid in pkgs
            FIX_PATCH[uuid] = true
            uuid in FIXED || push!(FIXED, uuid)
        end
    elseif opt == :fix_minor
        for uuid in pkgs
            FIX_MINOR[uuid] = true
            uuid in FIXED || push!(FIXED, uuid)
        end
    elseif opt == :fix_major
        for uuid in pkgs
            FIX_MAJOR[uuid] = true
            uuid in FIXED || push!(FIXED, uuid)
        end
    elseif opt == :unfix
        for uuid in pkgs,
            dict in (FIX_PATCH, FIX_MINOR, FIX_MAJOR)
            dict[uuid] = false
        end
        filter!(∉(pkgs), FIXED)
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

## options: sorting packages (resolution priority)

function default_sort_packages_by(uuid::UUID)
    uuid == JULIA_UUID ? ZERO_UUID : uuid
end

sort!(reqs, by = default_sort_packages_by)
sort!(FIXED, by = default_sort_packages_by)

sort_packages_by(uuid::UUID) = default_sort_packages_by(uuid)

let i = 0
    priority = Dict{UUID,Int}()
    handle_opts(:prioritize) do str::String
        for uuid in parse_packages(str)
            priority[uuid] = (i += 1)
        end
    end
    for uuid in FIXED
        priority[uuid] = (i += 1)
    end
    if !isempty(priority)
        n = i + 1
        global function sort_packages_by(uuid::UUID)
            get(priority, uuid, n), default_sort_packages_by(uuid)
        end
    end
end

## do an actual resolve

reg = registry_provider(
    packages;
    julia_versions,
    sort_versions,
    allow_pre,
    workspace_pkgs,
)

# the project's compat bounds are user constraints, not part of the package
# universe: they go into the `Problem`, which forbids the versions they exclude
# by clause instead of deleting them from the provider's data. (the `julia`
# bound is not among them: it was consumed above, as the Julia version
# universe.)
#
# bounds on the stdlibs Julia *pins* are widened to admit the bundled versions,
# so that answers stay exactly what they were when the provider applied compat
# itself: it applied compat to a package's registry versions and *then* patched
# the bundled ones back in, so those were never subject to it. a pinned stdlib's
# version is dictated by the chosen Julia, so a bound excluding it would only
# make that Julia infeasible -- Pkg treats such a bound as inert too.
#
# the upgradable stdlibs are excluded from the widening: nothing pins them, so
# their versions are resolved like any other package's, and a user bound on one
# must be enforced strictly -- again matching Pkg.
let compat = Dict{UUID,VersionSpec}(), bundled = bundled_versions(julia_versions, allow_pre)
    for (uuid, spec) in project_compat
        if uuid ∉ UPGRADABLE_STDLIBS_UUIDS
            for v in get(bundled, uuid, ())
                spec = spec ∪ VersionSpec(v)
            end
        end
        compat[uuid] = spec
    end
    global const problem = Problem(reqs; compat)
end

pkg_info = Resolver.pkg_info(reg, problem)
sol = resolve(pkg_info, problem; by=sort_packages_by)
sol === nothing && error("Unsatisfiable")

## output results

const julia_version = sol[JULIA_UUID]
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

# On Julia >= 1.14 the registry parser yields a package's dependencies as a
# `Set{UUID}` rather than a `Dict{String,UUID}`; the manifest still needs them
# keyed by name, so reconstruct the mapping from a uuid -> name table covering
# registered packages, every stdlib, and julia itself.
const uuid_to_name = Dict{UUID,String}(JULIA_UUID => "julia")
for (u, entries) in packages, entry in entries
    uuid_to_name[u] = entry.name
end
for (_v, stdlib_set) in STDLIBS_BY_VERSION, (u, info) in stdlib_set
    uuid_to_name[u] = info.name
end
for (u, info) in UNREGISTERED_STDLIBS
    uuid_to_name[u] = info.name
end
name_uuid_dict(d::AbstractDict) = d # <=1.13: already name => uuid
name_uuid_dict(d) = Dict{String,UUID}(uuid_to_name[u] => u for u in d) # >=1.14

const info_map = Dict{UUID,ManifestEntry}()

# is `uuid` recorded as a bundled stdlib of the resolved Julia, rather than as
# an ordinary registry package? every pinned stdlib is; an upgradable one only
# when it actually resolved to the version Julia bundles, since nothing pins it
# and a registry version of it is an ordinary package -- tree hash and registry
# deps -- exactly as Pkg records it. (`thispatch` drops the synthetic build
# numbers the provider adds to tell same-numbered stdlib entries apart.)
function bundled_stdlib(uuid::UUID, version::VersionNumber)
    uuid in keys(stdlibs) || return false
    uuid in UPGRADABLE_STDLIBS_UUIDS || return true
    bundled = something(stdlibs[uuid].version, julia_version)
    return Base.thispatch(version) == Base.thispatch(bundled)
end

for (uuid, version) in sol
    uuid === JULIA_UUID && continue
    # workspace packages become path entries in the manifest (written in the
    # manifest block below); they have no registry or stdlib info to record
    uuid in workspace_uuids && continue
    if bundled_stdlib(uuid, version)
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
        infos = Set{ManifestEntry}()
        for entry in packages[uuid]
            info = package_info(entry)
            haskey(info.version_info, version) || continue
            tree = info.version_info[version].git_tree_sha1
            deps = Dict{String,UUID}()
            for (r, d) in info.deps
                version in r || continue
                merge!(deps, name_uuid_dict(d))
            end
            if get(deps, "julia", nothing) == JULIA_UUID
                delete!(deps, "julia")
            end
            weakdeps = Dict{String,UUID}()
            for (r, d) in info.weak_deps
                version in r || continue
                merge!(weakdeps, name_uuid_dict(d))
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
    # the best version of a package that the project's constraints admit: the
    # info keeps the versions the Problem forbids (they are constrained away,
    # not deleted), so they don't count as available here
    function best_version(uuid::UUID)
        vers = pkg_info[uuid].versions
        i = findfirst(v -> !Resolver.is_excluded(problem, uuid, v), vers)
        return isnothing(i) ? first(vers) : vers[i]
    end
    # print packages and versions in priority order, required packages first
    pkgs = sort!(collect(keys(sol)), by = sort_packages_by)
    sort!(pkgs, by = !in(reqs))
    width = maximum(textwidth(info.name) for info in values(info_map))
    for uuid in pkgs
        uuid in workspace_uuids || continue
        global width = max(width, textwidth(workspace_pkgs[uuid][1]))
    end
    for uuid in pkgs
        if uuid == JULIA_UUID
            name = "julia"
            version = julia_version
        elseif uuid in workspace_uuids
            # workspace packages are fixed at their local version
            name, version, _ = workspace_pkgs[uuid]
        else
            name = info_map[uuid].name
            version = something(info_map[uuid].version, julia_version)
        end
        # a pinned stdlib had no choice to make; an upgradable one did
        optimal = uuid in keys(stdlibs) && uuid ∉ UPGRADABLE_STDLIBS_UUIDS ||
            version == best_version(uuid)
        try print(uuid, " ", rpad(name, width), " ", version)
            optimal || print(" ⊼")
            println()
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
    # add workspace projects as fixed path entries
    #
    # A workspace shares one manifest across the root project and its member
    # projects (Julia >= 1.12). Pkg records every workspace project that is a
    # *package* (i.e. has a uuid) as a `path` entry in that manifest -- the root
    # as `path = "."` and each member relative to the root -- so mirror that.
    # Bare sub-environments (test/docs/benchmarks with no uuid) contribute their
    # dependencies to resolution but are not themselves manifest entries.
    if in_workspace
        # `realpath` both ends: EnvCache stores workspace member paths resolved
        # of symlinks, so the base must be resolved the same way or `relpath`
        # would climb through the difference (e.g. /var vs /private/var on macOS,
        # or 8.3 short names on Windows) instead of yielding a clean relative path
        root_dir = realpath(dirname(PROJ))
        # the root project plus every member project, as (project-file, project)
        ws_projects = vcat([(PROJ, env.project)],
                           [(file, proj) for (file, proj) in env.workspace])
        # uuids allowed in a manifest entry's deps: everything resolved, plus the
        # workspace packages themselves (so inter-project deps are preserved)
        manifest_uuids = union(Set{UUID}(keys(deps)), workspace_uuids)
        for (proj_file, proj) in ws_projects
            proj.uuid === nothing && continue # bare sub-environment, not a package
            entry_deps = Dict{String,UUID}(
                name => uuid for (name, uuid) in proj.deps if uuid in manifest_uuids)
            entry_weakdeps = Dict{String,UUID}(
                name => uuid for (name, uuid) in proj.weakdeps if uuid in manifest_uuids)
            deps[proj.uuid] = PackageEntry(;
                uuid = proj.uuid,
                name = proj.name,
                version = proj.version,
                path = relpath(realpath(dirname(proj_file)), root_dir),
                deps = entry_deps,
                weakdeps = entry_weakdeps,
            )
        end
    end
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
