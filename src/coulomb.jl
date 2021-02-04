


using LinearAlgebra: norm

# this will contain all functions for implementing Ewald summation

function getChargeMatrix(charges::Array)
end


# Returns a list of lattice (real or reciprocal) vectors which are to be summed over. Equivalent of _buildList in coulomb.py
function getLatticeSummands(latticeVectors::Array, sumDepth::Int)
end


# This is the bulk ewald method
function ewald(q::Vector, Δ::Vector, crystal::Crystal, charges::Array)
end


# this is the slab ewald (deWette) method
function ewald(q::Vector, Δ::Vector, crystal::Slab, charges::Array)
end


function differentPlaneSum(q, Δ∥::Vector, Δ⟂::Vector, crystal::Slab)
end


# Equivalent of _DeWette in coulomb.py
function samePlaneSumDeWette(q::Vector, Δ::Vector, crystal::Slab)
end
