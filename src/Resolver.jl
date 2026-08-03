module Resolver

export resolve, issatisfiable,
    Problem, Diagnosis, Holdback, holdbacks

include("BitKernels.jl")
include("DepsProvider.jl")
include("Problem.jl")
include("PkgInfo.jl")
include("PkgInfoFiles.jl")
include("FilterPkgs.jl")
include("Goal.jl")
include("PicoSAT.jl")
include("SAT.jl")
include("Diagnose.jl")
include("Resolve.jl")
include("Closure.jl")
include("Holdback.jl")

end # module
