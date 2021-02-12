using LinearAlgebra: norm
# All the structs are in structure.jl

# this will contain all functions for implementing Ewald summation

function getChargeMatrix(charges::Array)
    "Get charge for matrix Z, returning a diagonal
    matrix containing the charges"
    e = 15.1891
    Z = Diagonal(repeat(charges, inner=[3])) .* e
end


# Returns a list of lattice (real or reciprocal) vectors which are to be summed over. Equivalent of _buildList in coulomb.py
function getLatticeSummands(latticeVectors::Array, sumDepth::Int)
    "Builds and returns a list of vectors to be summed over in Ewald summation"
    sumRange = collect(-sumDepth:1:sumDepth)   
    zSumRange = 
    if length(latticeVectors[1]) == 3
        
    
end


# This is the bulk ewald method
function ewald(q::Vector, Δ::Vector, crystal::Crystal, charges::Array)
end


# this is the slab ewald (deWette) method
function ewald(q::Vector, Δ::Vector, crystal::Slab, charges::Array)
end


function differentPlaneSum(q, Δparallel::Vector, Δnormal::Vector, crystal::Slab)
end


# Equivalent of _DeWette in coulomb.py
function samePlaneSumDeWette(q::Vector, Δ::Vector, crystal::Slab)
end
