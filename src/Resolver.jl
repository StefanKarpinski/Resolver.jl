module Resolver

export registry_provider, resolve
export PkgData, PkgInfo, load_pkg_info, make_pkg_info, filter_pkg_info!

include("DepsProvider.jl")
include("RegistryProvider.jl")
include("LoadPkgInfo.jl")
include("PkgInfoFiles.jl")
include("FilterPkgs.jl")
include("PicoSAT.jl")
include("Resolve.jl")

end # module
