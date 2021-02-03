

module LatticeDynamics

include("structure.jl")
include("shortrange.jl")
include("coulomb.jl")
include("dynamicalmatrix.jl")
include("model.jl")


export Crystal, Slab, getNeighbors!,
       buildPath, getDispersion, plotDispersion

end #module
