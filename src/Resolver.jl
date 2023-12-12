module Resolver

export registry_provider, resolve
export DepsProvider, PkgData, PkgInfo, SAT
export pkg_info, filter_pkg_info!
export
    finalize, sat_add, is_satisfiable,
    get_solution_inds!, with_temp_clauses, optimize_solution

include("DepsProvider.jl")
include("RegistryProvider.jl")
include("PkgInfo.jl")
include("PkgInfoFiles.jl")
include("FilterPkgs.jl")
include("PicoSAT.jl")
include("SAT.jl")
include("Resolve.jl")

end # module
