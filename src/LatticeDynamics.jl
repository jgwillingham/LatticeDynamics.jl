

module LatticeDynamics

using Distributed
using LinearAlgebra
using Printf
using GSL: sf_erfc, sf_gamma_inc
using ProgressMeter
using Plots


include("structure.jl")
include("shortrange.jl")
include("coulomb.jl")
include("dynamicalmatrix.jl")
include("greensfunction.jl")
include("model.jl")


export Crystal, Slab, projectVector,
       getSlabCouplingArray, getSlabCharges, buildPath,
       getDispersion, getProjectedDispersion,
       plotDispersion, 𝔻, φ,
       getSpectrum, plotSpectrum,
       getEnergySurface, plotEnergySurface

end #module
