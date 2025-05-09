module Registries

export registry_provider

import Base: UUID
import HistoricalStdlibVersions: STDLIBS_BY_VERSION, UNREGISTERED_STDLIBS, StdlibInfo
import JSON
import Pkg.Registry: JULIA_UUID, PkgEntry, init_package_info!, reachable_registries
import Pkg.Versions: VersionSpec
import Resolver: DepsProvider, PkgData

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

function registry_provider(;
    julia_versions :: VersionSpec = VersionSpec("1"),
    project_compat :: Dict{UUID,VersionSpec} = Dict{UUID,VersionSpec}(),
    sort_versions  :: Function = sort_versions_default,
    allow_pre      :: Dict{UUID,Bool} = Dict{UUID,Bool}(),
)
    packages = Dict{UUID,Vector{PkgEntry}}()

    for reg in reachable_registries()
        for (uuid, entry) in reg.pkgs
            push!(get!(()->PkgEntry[], packages, uuid), entry)
        end
    end

    function filter_pre!(uuid::UUID, vers::Vector{VersionNumber})
        if !get(allow_pre, uuid, allow_pre[UUID(0)])
            filter!(v->isempty(v.prerelease), vers)
        end
        return vers
    end

    julia_vers = filter(in(julia_versions), JULIA_VERSIONS)
    filter_pre!(JULIA_UUID, julia_vers)

    stdlibs = Dict{UUID, # stdlib UUID
        Dict{VersionNumber, # stdlib version
            Tuple{
                StdlibInfo, # stdlib info
                Vector{VersionNumber}, # compatible julia versions
            }
        }
    }()
    for julia_ver in julia_vers
        last_stdlibs = UNREGISTERED_STDLIBS
        for (v, this_stdlibs) in STDLIBS_BY_VERSION
            v ≥ Base.thispatch(julia_ver) && break
            last_stdlibs = this_stdlibs
        end
        for (uuid, stdlib_info) in last_stdlibs
            stdlib_ver = something(stdlib_info.version, julia_ver)
            stdlibs_u = get!(()->valtype(stdlibs)(), stdlibs, uuid)
            if stdlib_ver in keys(stdlibs_u)
                stdlibs_uv = stdlibs_u[stdlib_ver]
                @assert stdlibs_uv[1].name == stdlib_info.name
                @assert stdlibs_uv[1].uuid == stdlib_info.uuid
                if stdlibs_uv[1].version !== nothing
                    @assert stdlibs_uv[1].version == stdlib_ver
                end
                # @assert stdlibs_uv[1].deps == stdlib_info.deps
                if stdlibs_uv[1].deps != stdlib_info.deps
                    union!(stdlibs_uv[1].deps, stdlib_info.deps)
                end
                @assert stdlibs_uv[1].weakdeps == stdlib_info.weakdeps
                push!(stdlibs_uv[2], julia_ver)
            else
                stdlibs_u[stdlib_ver] = (stdlib_info, [julia_ver])
            end
        end
    end

    rp = DepsProvider(keys(packages)) do uuid::UUID
        vers = Set{VersionNumber}()
        deps = Dict{VersionNumber,Vector{UUID}}()
        comp = Dict{VersionNumber,Dict{UUID,VersionSpec}}()
        if uuid == JULIA_UUID
            union!(vers, julia_vers)
            for v in vers
                deps[v] = valtype(deps)()
                comp[v] = valtype(comp)()
            end
        elseif uuid in keys(packages)
            for entry in packages[uuid]
                info = init_package_info!(entry)
                # versions from this registry, filtered
                new_vers = collect(keys(info.version_info))
                filter_pre!(uuid, new_vers)
                if uuid in keys(project_compat)
                    filter!(in(project_compat[uuid]), new_vers)
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
                    # TODO: handle weakdeps
                end
            end
        elseif uuid ∉ keys(stdlibs)
            error("unknown package UUID: $uuid")
        end
        # patch in historical stdlib info
        if uuid in keys(stdlibs)
            for (v, (stdlib_info, julia_compat)) in stdlibs[uuid]
                v in vers && continue # prefer real registry data
                push!(vers, v)
                deps[v] = stdlib_info.deps
                comp[v] = Dict(JULIA_UUID => VersionSpec(julia_compat))
                # TODO: handle weakdeps
            end
        end
        # insert dependency on julia itself
        if uuid != JULIA_UUID
            for v in vers
                if JULIA_UUID ∉ deps[v]
                    push!(deps[v], JULIA_UUID)
                end
                if JULIA_UUID ∉ keys(comp[v])
                    comp[v][JULIA_UUID] = VersionSpec("*")
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

    return packages, stdlibs, rp
end

end # module

using .Registries
