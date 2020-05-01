

module LatticeDynamics

include("structure.jl")
include("model.jl")
include("dynamicalmatrix.jl")


export Crystal, getNeighbors!,
       buildPath, getDispersion, plotDispersion

end #module
