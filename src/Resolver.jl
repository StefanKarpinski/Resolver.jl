module Resolver

export registry_provider, resolve

include("DepsProvider.jl")
include("RegistryProvider.jl")
include("LoadPkgInfo.jl")
include("PkgInfoFiles.jl")
include("FilterPkgs.jl")
include("PicoSAT.jl")
include("Resolve.jl")

end # module
