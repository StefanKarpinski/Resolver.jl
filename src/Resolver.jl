module Resolver

export registry_provider, resolve
export DepsProvider, PkgData, PkgInfo, SAT, PicoSAT
export
    pkg_info, finalize, sat_add, sat_assume, is_satisfiable, is_unsatisfiable,
    extract_solution!, optimize_solution!, with_temp_clauses

include("DepsProvider.jl")
include("RegistryProvider.jl")
include("PkgInfo.jl")
include("PkgInfoFiles.jl")
include("FilterPkgs.jl")
include("PicoSAT.jl")
include("SAT.jl")
include("Resolve.jl")

end # module
