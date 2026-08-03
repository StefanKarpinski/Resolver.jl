module Resolver

export resolve, Problem, Diagnosis

include("BitKernels.jl")
include("DepsProvider.jl")
include("Problem.jl")
include("PkgInfo.jl")
include("PkgInfoFiles.jl")
include("FilterPkgs.jl")
include("PicoSAT.jl")
include("SAT.jl")
include("Diagnose.jl")
include("Resolve.jl")
include("Closure.jl")

end # module
