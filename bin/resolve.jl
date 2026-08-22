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
  --allow-yanked[=<pkgs>] allow yanked versions
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
    julia allow-pre allow-yanked extra-deps prioritize
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
import Resolver.Diagnostics
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
# Whichever source wins is an ordinary compat entry on `julia`: the provider
# offers every Julia version there is, along with every stdlib version any of
# them bundles (see Registries.jl), and this bound constrains the `julia`
# package the way any other bound constrains its own. So a project bound and
# the equivalent `--julia` are indistinguishable downstream -- and the Julias
# the bound excludes are excluded by clause rather than missing from the
# universe, which is what lets one artifact answer for any of them. The
# stdlib <-> Julia couplings do the rest: a Julia the bound rules out cannot
# be chosen, so neither can a stdlib version only that Julia bundles.

const julia_versions = something(
    handle_opts(:julia) do val::String
        try VersionSpec(split(val, r"\s*,\s*"))
        catch
            usage("Invalid compat version spec: --julia=$val")
        end
    end,
    get(project_compat, JULIA_UUID, nothing),
    VersionSpec("1"),
)

project_compat[JULIA_UUID] = julia_versions

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
const allow_yanked = Dict(ZERO_UUID => false)
const FIX_PATCH = Dict(ZERO_UUID => false)
const FIX_MINOR = Dict(ZERO_UUID => false)
const FIX_MAJOR = Dict(ZERO_UUID => false)
const ORDER_MAP = Dict(ZERO_UUID => :max)

# list of all fixed packages
const FIXED = UUID[]

handle_opts(r"^(allow_(pre|yanked)|(un)?fix|max|min)") do opt, val
    pkgs = isnothing(val) ? [ZERO_UUID] : parse_packages(val)
    if opt == :allow_pre
        for uuid in pkgs
            allow_pre[uuid] = true
        end
    elseif opt == :allow_yanked
        for uuid in pkgs
            allow_yanked[uuid] = true
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

## version ordering
#
# The version preference ordering is a *query* parameter, not part of the package
# universe: it decides which of the valid solutions is optimal, not which
# solutions are valid. So the provider hands over versions in the canonical order
# (newest first, see Registries.jl) and this ordering reaches the resolver as
# `resolve`'s `order` argument -- one comparator per package, built from the
# --fix*/--max*/--min* options and the old manifest.

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

function level_order(level::Symbol) :: Function
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

function version_order(uuid::UUID) :: Function
    fixed_by = fixed(
        get(old_versions, uuid, nothing),
        get(FIX_PATCH, uuid, FIX_PATCH[ZERO_UUID]),
        get(FIX_MINOR, uuid, FIX_MINOR[ZERO_UUID]),
        get(FIX_MAJOR, uuid, FIX_MAJOR[ZERO_UUID]),
    )
    <ₒ = level_order(get(ORDER_MAP, uuid, ORDER_MAP[ZERO_UUID]))
    function lt(u::VersionNumber, v::VersionNumber)
        fixed_u = fixed_by(u) :: UInt8
        fixed_v = fixed_by(v) :: UInt8
        fixed_u ≠ fixed_v ? fixed_u > fixed_v : u <ₒ v
    end
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

# The versions the yanked filter dropped, per package, filled in by the provider
# as it is asked for packages. Yankedness is registry-derived, not a query knob:
# the universe is built without the versions the registries have struck, so a
# resolve can neither produce one nor suggest one. `--allow-yanked` asks for them
# back -- for the named packages, or for all of them when given bare -- which
# builds a universe that keeps them.
const yanked_dropped = Dict{UUID,Set{VersionNumber}}()

reg = registry_provider(packages;
    workspace_pkgs, allow_yanked, yanked = yanked_dropped)

# the project's compat bounds are user constraints, not part of the package
# universe: they go into the `Problem`, which forbids the versions they exclude
# by clause instead of deleting them from the provider's data. the `julia` bound
# is one of them, no longer special (see above).
#
# every bound is enforced strictly, stdlibs included. a bound on a stdlib is not
# inert: a Julia version is compatible with exactly the stdlib versions it
# bundles -- that is what its compat pins say -- so a bound that excludes those
# versions excludes the Julia along with them, and the resolver *steers* the
# Julia choice to satisfy it. If no admissible Julia bundles a version the bound
# admits, and no registry version fits either, the requirements really are
# unsatisfiable and saying so beats silently ignoring the bound.
#
# prerelease admission rides along as an exclusion kind: `--allow-pre` says which
# packages' prereleases the query accepts, the same way a compat bound forbids a
# version. (Yankedness is not one of these: it is a property of the registry state,
# so the struck versions are not in the universe to begin with -- see the provider.)
const problem = Problem(reqs;
    compat = project_compat,
    prerelease = prerelease_exclusion(allow_pre),
)

pkg_info = Resolver.pkg_info(reg, problem)
sol = resolve(pkg_info, problem; by=sort_packages_by, order=version_order)

# The one unsatisfiable shape the universe cannot explain for itself: a bound
# whose only admissible versions are ones the yanked filter dropped. Those
# versions are plainly there in the registry, so "unsatisfiable" on its own reads
# as a lie -- name them, and name the flag that brings them back.
function yanked_notes()
    notes = String[]
    for uuid in sort!(collect(keys(project_compat)))
        uuid in keys(packages) || continue
        spec = project_compat[uuid]
        # asking the provider for the package both fills in what it dropped and
        # says what it kept, so the note doesn't depend on whether the resolve
        # happened to reach this package
        kept = Resolver.pkg_data(reg, uuid).versions
        dropped = get(yanked_dropped, uuid, nothing)
        dropped === nothing && continue
        admitted = sort!([v for v in dropped if v ∈ spec])
        isempty(admitted) && continue
        # only when nothing that survived the filter is admissible either:
        # otherwise the bound is satisfiable and the yank is not the problem
        any(v -> v ∈ spec, kept) && continue
        name = first(packages[uuid]).name
        push!(notes, "note: the versions your compat on $name admits are " *
            "yanked ($(join(admitted, ", "))); --allow-yanked accepts them anyway")
    end
    return notes
end

# A diagnosis is plain data, keyed by whatever the universe is keyed by -- here
# uuids, which is not what the reader knows the packages as. So rebuild it over
# names before printing it: exactly the "render your own report" the type is
# shaped for, minus the rendering, since the default one is what we want.

function named(d::Resolver.Diagnosis, name::Function)
    Resolver.Diagnosis(
        [Diagnostics.Conflict{String,VersionNumber}(
            String[name(p) for p in c.reqs],
            Diagnostics.Fact[named(f, name) for f in c.chain],
            [Diagnostics.Fix{String,VersionNumber}(
                [Diagnostics.Action(a.kind, name(a.pkg)) for a in fix.actions],
                Dict{String,VersionNumber}(
                    name(p) => v for (p, v) in fix.solution))
             for fix in c.fixes])
         for c in d.conflicts],
        d.others)
end

named(f::Diagnostics.Requirement, name::Function) =
    Diagnostics.Requirement(name(f.pkg))
named(f::Diagnostics.Availability, name::Function) =
    Diagnostics.Availability(name(f.pkg), f.members, f.excluded)
named(f::Diagnostics.Dependency, name::Function) =
    Diagnostics.Dependency(name(f.pkg), f.versions, name(f.dep),
        f.offering, f.allowed)
named(f::Diagnostics.Incompatibility, name::Function) =
    Diagnostics.Incompatibility(name(f.pkg), f.versions, name(f.other),
        f.offering, f.allowed)

function package_name(uuid::UUID)
    uuid == JULIA_UUID && return "julia"
    haskey(workspace_pkgs, uuid) && return workspace_pkgs[uuid][1]
    haskey(packages, uuid) && return first(packages[uuid]).name
    for (_, this_stdlibs) in STDLIBS_BY_VERSION
        haskey(this_stdlibs, uuid) && return this_stdlibs[uuid].name
    end
    haskey(UNREGISTERED_STDLIBS, uuid) && return UNREGISTERED_STDLIBS[uuid].name
    return string(uuid)
end

if sol isa Resolver.Diagnosis
    show(stderr, MIME("text/plain"), named(sol, package_name))
    for note in yanked_notes()
        println(stderr, note)
    end
    exit(1)
end

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
    # The best version of a package that was actually available. The universe
    # holds every version of every package and every Julia, and says what a query
    # admits by constraining rather than by leaving things out -- so "available"
    # has to be spelled out here: the project's bounds must admit it, and it must
    # be a version the registries have or one that a Julia the bound admits
    # bundles. (A version bundled only by some other Julia is in the universe but
    # was never on offer, and would make every stdlib look sub-optimal.) The info
    # lists versions in the canonical order, not the requested one, so rank them
    # with the same comparator the resolve used.
    bundled_here = bundled_versions(julia_versions)
    function on_offer(uuid::UUID, v::VersionNumber)
        # a version the registries have is on offer whatever Julia we target,
        # and so is anything no Julia bundles (`julia` itself included)
        is_bundled(uuid, v) || return true
        is_registered(packages, uuid, v) && return true
        # otherwise it exists only because some Julia bundles it
        return v in get(bundled_here, uuid, ())
    end
    function best_version(uuid::UUID)
        vers = sort!(copy(pkg_info[uuid].versions); lt = version_order(uuid))
        i = findfirst(v -> !Resolver.is_excluded(problem, uuid, v) &&
                           on_offer(uuid, v), vers)
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
        ctx = Context(; env, julia_version)
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
