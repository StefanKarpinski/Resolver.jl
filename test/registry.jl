using Pkg: depots1
using Pkg.Registry: PkgInfo, RegistryInstance, init_package_info!
using Pkg.Types: stdlibs
using Pkg.Versions: VersionSpec

const EXCLUDED_DEPS = push!(Set(values(stdlibs())), "julia")

const reg_path = "~/.julia/registries/General.toml"
const reg = RegistryInstance(expanduser(reg_path))

const PKG_INFO = Dict{String,PkgInfo}()

pkg_info(p::AbstractString) = get!(PKG_INFO, p) do
    for pkg in values(reg.pkgs)
        pkg.name == p && return init_package_info!(pkg)
    end
    error("package $p not found")
end

function versions_callback(p::AbstractString)
    info = pkg_info(p)
    sort!(collect(keys(info.version_info)), rev=true)
end

function deps_callback(p::AbstractString, v::VersionNumber)
    info = pkg_info(p)
    deps = Set{String}()
    for (r, d) in info.deps
        v in r && union!(deps, keys(d))
    end
    return setdiff(deps, EXCLUDED_DEPS)
end

function compat_callback(p₁::AbstractString, v₁::VersionNumber)
    info = pkg_info(p₁)
    compat = Dict{String,VersionSpec}()
    for (r₁, c₁) in info.compat
        v₁ in r₁ && merge!(compat, c₁)
    end
    function (p₂::AbstractString, v₂::VersionNumber)
        r₂ = get(compat, p₂, nothing)
        r₂ === nothing || v₂ in r₂
    end
end
