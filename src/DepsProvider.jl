# Type Parameters:
#  P = package type (eg String)
#  V = version type (eg VersionNumber)
#  S = version set type (eg VersionSpec)

const SetOrVec{T} = Union{AbstractSet{T}, AbstractVector{T}}

struct PkgData{P,V,S}
    versions :: Vector{V}
    depends  :: Dict{V, Vector{P}}
    compat   :: Dict{V, Dict{P, S}}
end

struct DepsProvider{P,V,S,F<:Function}
    packages :: Vector{P}
    provider :: F
end

function DepsProvider{P,V,S}(
    provider :: Function,
    packages :: SetOrVec{P},
) where {P,V,S}
    packages = sort!(P[p for p in packages])
    DepsProvider{P,V,S,typeof(provider)}(packages, provider)
end

(deps::DepsProvider{P,V,S,F})(pkg::P) where {P,V,S,F<:Function} =
    deps.provider(pkg) :: PkgData{P,V,S}
