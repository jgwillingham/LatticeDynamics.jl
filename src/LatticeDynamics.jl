

module LatticeDynamics

include("structure.jl")
include("shortrange.jl")
include("coulomb.jl")
include("dynamicalmatrix.jl")
include("greensfunction.jl")
include("model.jl")


export Crystal, Slab, projectVector,
       getSlabCouplingArray, buildPath,
       getDispersion, getProjectedDispersion,
       plotDispersion, 𝔻,
       getSpectrum, plotSpectrum,
       getEnergySurface, plotEnergySurface

end #module
