module Registries

export registry_provider, package_info, bundled_versions,
    UPGRADABLE_STDLIBS_UUIDS

import Base: UUID
import HistoricalStdlibVersions: STDLIBS_BY_VERSION, UNREGISTERED_STDLIBS, StdlibInfo
import JSON
import Pkg
import Pkg.Registry: JULIA_UUID, PkgEntry, RegistryInstance, init_package_info!, isyanked, reachable_registries
import Pkg.Versions: VersionSpec
import Resolver: DepsProvider, PkgData

## the upgradable stdlibs
#
# Julia bundles a version of every stdlib, but two of them it does not *pin*:
# `Pkg.Types.UPGRADABLE_STDLIBS`. For those, registry versions compete with the
# bundled one like any other package's, which is what Pkg does and what we
# model (see the JULIA_UUID branch of the provider). Pkg only grew the constant
# in 1.12, so fall back to its values on the older Pkgs `bin/` supports.
const UPGRADABLE_STDLIBS_UUIDS =
    if isdefined(Pkg.Types, :UPGRADABLE_STDLIBS_UUIDS)
        Set{UUID}(Pkg.Types.UPGRADABLE_STDLIBS_UUIDS)
    else
        Set{UUID}([
            UUID("8bb1440f-4735-579b-a4ab-409b98df4dab"), # DelimitedFiles
            UUID("10745b16-79ce-11e8-11f9-7d13ad32a3b2"), # Statistics
        ])
    end

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

# share one object among equal-valued entries; there are only a few
# distinct values per package, so scanning them beats comparing all pairs
function dedup_values!(d::AbstractDict, vers::Vector)
    distinct = valtype(d)[]
    for v in vers
        x = d[v]
        shared = false
        for y in distinct
            if y == x
                d[v] = y
                shared = true
                break
            end
        end
        shared || push!(distinct, x)
    end
    return d
end

# The Julia versions to resolve against, and what they bundle:
#
#   stdlibs[uuid][version] : the historical stdlib info for a bundled version
#   stdlib_pins[uuid]      : the (Julia version, bundled-version spec) pairs.
#     Recorded so `registry_provider` can reconcile a stdlib version's `julia`
#     compat with the Julia versions that bundle it. Recorded for the
#     upgradable stdlibs too, which are bundled without being pinned: a bundled
#     version must be installable on its Julia either way.
#
# Split out of `registry_provider` so that `bundled_versions` can report the
# bundled version sets without building a provider.
function julia_and_stdlib_versions(
    julia_versions :: VersionSpec,
    allow_pre      :: Dict{UUID,Bool} = Dict(UUID(0) => false),
)
    julia_vers = filter(in(julia_versions), JULIA_VERSIONS)
    get(allow_pre, JULIA_UUID, allow_pre[UUID(0)]) ||
        filter!(v -> isempty(v.prerelease), julia_vers)

    stdlibs = Dict{UUID,Dict{VersionNumber,StdlibInfo}}()
    stdlib_pins = Dict{UUID,Vector{Tuple{VersionNumber,VersionSpec}}}()
    for julia_ver in julia_vers
        last_stdlibs = UNREGISTERED_STDLIBS
        for (v, this_stdlibs) in STDLIBS_BY_VERSION
            v > Base.thispatch(julia_ver) && break
            last_stdlibs = this_stdlibs
        end
        for (uuid, stdlib_info) in last_stdlibs
            stdlib_ver = something(stdlib_info.version, julia_ver)
            push!(get!(()->valtype(stdlib_pins)(), stdlib_pins, uuid),
                  (julia_ver, VersionSpec(stdlib_ver)))
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
    return julia_vers, stdlibs, stdlib_pins
end

# Per package, the versions the provider offers because a Julia in
# `julia_versions` bundles them, regardless of what the registries say. Nothing
# in the resolve path needs this -- a bundled version is an ordinary candidate,
# and a user bound that excludes it excludes the Julias that ship it -- but it
# is the answer to "which versions of this stdlib can I even get on these
# Julias", so it stays available for introspection and for the tests.
function bundled_versions(
    julia_versions :: VersionSpec,
    allow_pre      :: Dict{UUID,Bool} = Dict(UUID(0) => false),
)
    _, stdlibs, _ = julia_and_stdlib_versions(julia_versions, allow_pre)
    Dict{UUID,Vector{VersionNumber}}(
        uuid => collect(keys(vers)) for (uuid, vers) in stdlibs)
end

function registry_provider(
    packages       :: Dict{UUID,Vector{PkgEntry}};
    julia_versions :: VersionSpec = VersionSpec("1"),
    sort_versions  :: Function = sort_versions_default,
    allow_pre      :: Dict{UUID,Bool} = Dict{UUID,Bool}(),
    workspace_pkgs :: Dict{UUID,Tuple{String,VersionNumber,Vector{UUID}}} =
                      Dict{UUID,Tuple{String,VersionNumber,Vector{UUID}}}(),
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

    julia_vers, stdlibs, stdlib_pins =
        julia_and_stdlib_versions(julia_versions, allow_pre)

    return DepsProvider(keys(packages)) do uuid::UUID
        if uuid in keys(workspace_pkgs)
            # a workspace package is fixed at its local version: exactly one
            # version, with the local project's deps plus the implicit
            # dependency on julia itself. The registry is never consulted --
            # even when the uuid shadows a registered package the local copy
            # wins, and dependents' compat is enforced against the fixed
            # version by ordinary compat propagation, just as Pkg treats
            # fixed path dependencies.
            _, v, ws_deps = workspace_pkgs[uuid]
            deps_v = copy(ws_deps)
            JULIA_UUID in deps_v || push!(deps_v, JULIA_UUID)
            comp_v = Dict{UUID,VersionSpec}(JULIA_UUID => VersionSpec("*"))
            return PkgData([v], Dict(v => deps_v), Dict(v => comp_v))
        end
        vers = Set{VersionNumber}()
        deps = Dict{VersionNumber,Vector{UUID}}()
        comp = Dict{VersionNumber,Dict{UUID,VersionSpec}}()
        # versions that exist only because a Julia bundles them (no registry
        # entry); the widening below bounds them to their bundling Julias
        bundled_only = Set{VersionNumber}()
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
                # pin every stdlib this Julia bundles to its bundled version
                # -- except the upgradable ones, which Julia bundles without
                # pinning: for those the registry versions compete with the
                # bundled one on their own `julia` compat, exactly as Pkg
                # resolves them. "Pinned" here means pinned per candidate
                # Julia: we resolve over Julia versions too, so each Julia
                # version carries its own pins, and an upgradable stdlib
                # simply carries none.
                comp_v = get!(()->valtype(comp)(), comp, v)
                for (stdlib_uuid, stdlib_info) in last_stdlibs
                    stdlib_uuid in UPGRADABLE_STDLIBS_UUIDS && continue
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
                # the project's own compat is *not* applied here: it is a user
                # constraint, and it reaches the resolver as part of the
                # `Problem` (see bin/resolve.jl). the provider only decides
                # which versions exist at all
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
                push!(bundled_only, v)
                deps[v] = info.deps
            end
        end
        # widen a stdlib package's `julia` compat to admit the Julias that bundle it
        #
        # For a package that is also a stdlib, the Julia <-> stdlib version
        # coupling is expressed authoritatively by Julia's own compat pins (see
        # the JULIA_UUID branch above, where each Julia version pins every one
        # of its stdlibs to an exact version). A registry entry's own `julia`
        # compat for such a package can *disagree* with a Julia that actually
        # bundles a given version -- e.g. CompilerSupportLibraries_jll 1.1.1
        # ships with Julia 1.10.8+, yet the General registry marks 1.1.x as
        # requiring `julia = "1.11.0-1"`. Left as-is, that bound makes the
        # pinned stdlib version conflict with the very Julia that pins it, so
        # anything pulling in the BLAS stack (LinearAlgebra -> OpenBLAS_jll ->
        # CompilerSupportLibraries_jll) is unsatisfiable there.
        #
        # For a version that *is* in the registry we widen rather than replace:
        # it keeps its registry `julia` compat and additionally admits every
        # Julia that bundles it. This is important because the same version can
        # be a bundled stdlib for one Julia yet a resolvable dependency for
        # another; for the latter the registry bound must still govern, so
        # dropping it outright would wrongly make the version installable on
        # Julias it isn't compatible with.
        #
        # A version that exists *only* as a bundled stdlib has no registry bound
        # to keep, and it is not installable anywhere except where it ships:
        # its `julia` compat is exactly the Julias that bundle it. Leaving it
        # unbounded would let it pair with any Julia -- which the pins used to
        # rule out for pinned stdlibs, but nothing rules out for the upgradable
        # ones (so a Julia 1.12 resolve could take the Statistics that only
        # Julia 1.10 ships).
        #
        # The upgradable stdlibs get the same treatment even though no pin makes
        # them conflict: a bundled version must be installable on the Julia that
        # bundles it whether or not that Julia insists on it.
        if uuid in keys(stdlib_pins)
            # per version, the Julias that bundle it
            bundling = Dict{VersionNumber,VersionSpec}()
            for (julia_ver, ver_spec) in stdlib_pins[uuid]
                for v in vers
                    v in ver_spec || continue # julia_ver bundles this version
                    bundling[v] = haskey(bundling, v) ?
                        bundling[v] ∪ VersionSpec(julia_ver) : VersionSpec(julia_ver)
                end
            end
            for (v, julias) in bundling
                comp_v = get!(()->valtype(comp)(), comp, v)
                comp_v[JULIA_UUID] = v in bundled_only ? julias :
                    get(comp_v, JULIA_UUID, VersionSpec("*")) ∪ julias
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
        dedup_values!(deps, vers)
        dedup_values!(comp, vers)
        # return resolver PkgData structure
        PkgData(vers, deps, comp)
    end
end

end # module

using .Registries
