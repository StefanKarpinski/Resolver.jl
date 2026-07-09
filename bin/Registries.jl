module Registries

export registry_provider, package_info

import Base: UUID
import HistoricalStdlibVersions: STDLIBS_BY_VERSION, UNREGISTERED_STDLIBS, StdlibInfo
import JSON
import Pkg.Registry: JULIA_UUID, PkgEntry, RegistryInstance, init_package_info!, isyanked, reachable_registries
import Pkg.Versions: VersionSpec
import Resolver: DepsProvider, PkgData

## metaprogram around Pkg's `init_package_info!` signature change
#
# In Julia 1.14 the one-argument `init_package_info!(::PkgEntry)` was removed in
# favor of `init_package_info!(::RegistryInstance, ::PkgEntry)`. A PkgEntry only
# records its `registry_path`, so we map that back to the RegistryInstance that
# owns it. `package_info(entry)` works on both old and new Pkg.
if hasmethod(init_package_info!, Tuple{RegistryInstance, PkgEntry})
    const REGISTRIES_BY_PATH =
        Dict(reg.path => reg for reg in reachable_registries())
    package_info(entry::PkgEntry) =
        init_package_info!(REGISTRIES_BY_PATH[entry.registry_path], entry)
else
    package_info(entry::PkgEntry) = init_package_info!(entry)
end

## metaprogram around the Pkg 1.14 PkgInfo representation change
#
# Julia 1.14 changed the parsed registry data from name-keyed to uuid-keyed:
#   deps   value: Dict{String,UUID} (<=1.13)  ->  Set{UUID}            (>=1.14)
#   compat value: Dict{String,VersionSpec}     ->  Dict{UUID,VersionSpec}
# These helpers read either representation and yield uuid-based data, which is
# what we build regardless.

# the dependency uuids in a deps-map value
dep_uuids(d::AbstractDict) = values(d)   # <=1.13: name => uuid
dep_uuids(d::AbstractSet)  = d            # >=1.14: Set{UUID}

# (uuid => spec) pairs from a compat-map value; `name2uuid` resolves names on
# the old (name-keyed) representation and is ignored on the new one
compat_uuid_pairs(c::AbstractDict{<:AbstractString}, name2uuid::AbstractDict) =
    (name2uuid[n] => s for (n, s) in c)
compat_uuid_pairs(c::AbstractDict{UUID}, name2uuid::AbstractDict) = c

## download Julia versions

import Downloads

julia_versions_url = "https://julialang-s3.julialang.org/bin/versions.json"
julia_versions_file = joinpath(@__DIR__, "julia_versions.json")

if !isfile(julia_versions_file) ||
    time() - mtime(julia_versions_file) > 3600 # 1 hour
    Downloads.download(julia_versions_url, julia_versions_file)
end

julia_versions_data = JSON.parsefile(julia_versions_file)
const JULIA_VERSIONS = sort!(VersionNumber.(keys(julia_versions_data)))

# drop prerelease versions unless there is no corresponding final release
# when that happens only keep the last prerelease
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

function registry_provider(
    packages       :: Dict{UUID,Vector{PkgEntry}};
    julia_versions :: VersionSpec = VersionSpec("1"),
    project_compat :: Dict{UUID,VersionSpec} = Dict{UUID,VersionSpec}(),
    sort_versions  :: Function = sort_versions_default,
    allow_pre      :: Dict{UUID,Bool} = Dict{UUID,Bool}(),
)
    function filter_pre!(uuid::UUID, vers::Vector{VersionNumber})
        if !get(allow_pre, uuid, allow_pre[UUID(0)])
            filter!(v->isempty(v.prerelease), vers)
        end
        return vers
    end

    function filter_yanked!(info, vers::Vector{VersionNumber})
        filter!(v -> !isyanked(info, v), vers)
        return vers
    end

    julia_vers = filter(in(julia_versions), JULIA_VERSIONS)
    filter_pre!(JULIA_UUID, julia_vers)

    stdlibs = Dict{UUID,Dict{VersionNumber,StdlibInfo}}()
    for julia_ver in julia_vers
        last_stdlibs = UNREGISTERED_STDLIBS
        for (v, this_stdlibs) in STDLIBS_BY_VERSION
            v > Base.thispatch(julia_ver) && break
            last_stdlibs = this_stdlibs
        end
        for (uuid, stdlib_info) in last_stdlibs
            stdlib_ver = something(stdlib_info.version, julia_ver)
            deps_u = get!(()->valtype(stdlibs)(), stdlibs, uuid)
            # Sometimes HistoricalStdlibVersions gives the same
            # version number to stdlib entries with different deps.
            # This is horrible and bad and wrong. To deal with this
            # we tack on a build number to make different versions.
            # This is a horrible hack, but desperate times and all.
            build = 0
            original_ver = stdlib_ver
            while stdlib_ver in keys(deps_u) && !(
                    deps_u[stdlib_ver].name     == stdlib_info.name &&
                    deps_u[stdlib_ver].uuid     == stdlib_info.uuid &&
                    deps_u[stdlib_ver].deps     == stdlib_info.deps &&
                    deps_u[stdlib_ver].weakdeps == stdlib_info.weakdeps
                )
                stdlib_ver = VersionNumber(
                    original_ver.major,
                    original_ver.minor,
                    original_ver.patch,
                    original_ver.prerelease,
                    (original_ver.build..., "fake", build += 1),
                )
            end
            deps_u[stdlib_ver] = stdlib_info
        end
    end

    return DepsProvider(keys(packages)) do uuid::UUID
        vers = Set{VersionNumber}()
        deps = Dict{VersionNumber,Vector{UUID}}()
        comp = Dict{VersionNumber,Dict{UUID,VersionSpec}}()
        if uuid == JULIA_UUID
            union!(vers, julia_vers)
            for v in vers
                deps[v] = valtype(deps)()
                # find relevant stdlibs stanza
                last_stdlibs = UNREGISTERED_STDLIBS
                for (v′, this_stdlibs) in STDLIBS_BY_VERSION
                    v′ > Base.thispatch(v) && break
                    last_stdlibs = this_stdlibs
                end
                # add compat for all stdlibs of this version
                comp_v = get!(()->valtype(comp)(), comp, v)
                for (stdlib_uuid, stdlib_info) in last_stdlibs
                    stdlib_ver = something(stdlib_info.version, v)
                    comp_v[stdlib_uuid] = VersionSpec(stdlib_ver)
                end
            end
        elseif uuid in keys(packages)
            for entry in packages[uuid]
                info = package_info(entry)
                # versions from this registry, filtered
                new_vers = collect(keys(info.version_info))
                filter_pre!(uuid, new_vers)
                filter_yanked!(info, new_vers)
                if uuid in keys(project_compat)
                    filter!(in(project_compat[uuid]), new_vers)
                end
                # scan versions and populate deps & compat data
                for v in new_vers
                    # NOTE: we probably won't support the same name meaning
                    # different things in deps vs weakdeps, but here we can
                    # allow it, so we defensively do allow it
                    push!(vers, v)
                    # strong deps
                    deps_uuids = Dict{String,UUID}() # name => uuid (<=1.13 only)
                    deps_v = get!(()->valtype(deps)(), deps, v)
                    for (r, d) in info.deps
                        v in r || continue
                        d isa AbstractDict && merge!(deps_uuids, d)
                        union!(deps_v, dep_uuids(d))
                    end
                    # strong compat
                    comp_v = get!(()->valtype(comp)(), comp, v)
                    for (r, c) in info.compat
                        v in r || continue
                        for (u, spec) in compat_uuid_pairs(c, deps_uuids)
                            if u in keys(comp_v)
                                comp_v[u] = spec ∩ comp_v[u]
                            else
                                comp_v[u] = spec
                            end
                        end
                    end
                    # weak deps
                    weak_uuids = Dict{String,UUID}() # name => uuid (<=1.13 only)
                    for (r, d) in info.weak_deps
                        v in r || continue
                        d isa AbstractDict && merge!(weak_uuids, d)
                    end
                    # weak compat
                    for (r, c) in info.weak_compat
                        v in r || continue
                        for (u, spec) in compat_uuid_pairs(c, weak_uuids)
                            if u in keys(comp_v)
                                comp_v[u] = spec ∩ comp_v[u]
                            else
                                comp_v[u] = spec
                            end
                        end
                    end
                end
            end
        elseif uuid ∉ keys(stdlibs)
            error("unknown package UUID: $uuid")
        end
        # patch in historical stdlib info
        if uuid in keys(stdlibs)
            for (v, info) in stdlibs[uuid]
                v in vers && continue # prefer real registry data
                push!(vers, v)
                deps[v] = info.deps
            end
        end
        # insert dependency on julia itself
        if uuid != JULIA_UUID
            for v in vers
                deps_v = get!(()->valtype(deps)(), deps, v)
                if JULIA_UUID ∉ deps_v
                    push!(deps_v, JULIA_UUID)
                end
                comp_v = get!(()->valtype(comp)(), comp, v)
                if JULIA_UUID ∉ keys(comp_v)
                    comp_v[JULIA_UUID] = VersionSpec("*")
                end
            end
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
end

end # module

using .Registries
