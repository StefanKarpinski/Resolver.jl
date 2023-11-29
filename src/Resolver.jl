module Resolver

export
    registry_provider,
    load_pkg_info,
    save_pkg_info_file,
    load_pkg_info_file,
    find_reachable!,
    find_redundant!,
    shrink_pkg_info!,
    resolve

include("DepsProvider.jl")
include("RegistryProvider.jl")
include("LoadPkgInfo.jl")
include("SavePkgInfo.jl")
include("FilterPkgs.jl")
include("PicoSAT.jl")
include("Resolve.jl")

end # module
