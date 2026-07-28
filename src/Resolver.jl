module Resolver

export resolve

include("BitKernels.jl")
include("DepsProvider.jl")
include("PkgInfo.jl")
include("PkgInfoFiles.jl")
include("FilterPkgs.jl")
include("PicoSAT.jl")
include("SAT.jl")
include("Resolve.jl")

end # module
