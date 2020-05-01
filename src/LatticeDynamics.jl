

module LatticeDynamics

include("structure.jl")
include("shortrange.jl")
include("ewald.jl")
include("dynamicalmatrix.jl")
include("model.jl")


export Crystal, getNeighbors!,
       buildPath, getDispersion, plotDispersion

end #module
