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
    sumRange = -sumDepth:1:sumDepth
    zSumRange = Vector{}
    v = latticeVectors
    if length(latticeVectors[1]) == 3
        zSumRange = sumRange
    else if length(latticeVectors[1]) == 2
        zSumRange = [0]
        push!(v,zeros(3))
    else
        print("latticeVectors size inappropriate")
    end

    vec(n1,n2,n3) = n1*v[1] + n2*v[2] + n3*v[3]

    l = []
    for n1 in sumRange
        for n2 in sumRange
            for n3 in sumRange
                if n1==n2==n3
                    continue
                else
                    V = vec(n1,n2,n3)
                    append!(l, V)
                end
            end
        end
    end
    return l
end


# This is the bulk ewald method
function bulkEwald(q::Vector, Δ::Vector, crystal::Crystal, charges::Array)
    """
    Calculates the Ewald summation for the bulk crystal at wavevector `q`. 
    It returns the block of the Coulomb contribution to the dynamical
    matrix relating the two atoms separated by `INTRACELL_DISTANCE`

    Returns
    -------
    C_ij : matrix
        Block i,j of Coulomb contribution to dynamical matrix.

    """
    C_far = QSPACESUM(q, INTRACELL_DISTANCE)
    C_near = REALSPACESUM(q,INTRACELL_DISTANCE)
    C_ij = C_far + C_near
end


# this is the slab ewald (deWette) method
function slabEwald(q::Vector, Δ::Vector, crystal::Slab, charges::Array)
end


function differentPlaneSum(q, Δparallel::Vector, Δnormal::Vector, crystal::Slab)
end


# Equivalent of _DeWette in coulomb.py
function samePlaneSumDeWette(q::Vector, Δ::Vector, crystal::Slab)
end
