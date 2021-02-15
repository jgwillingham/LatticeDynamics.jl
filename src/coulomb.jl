using LinearAlgebra: norm

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
    v = latticeVectors
    if length(latticeVectors) == 3
        zSumRange = sumRange
    else if length(latticeVectors) == 2
        zSumRange = [0]
        push!(v,zeros(3))
    else
        print("latticeVectors size inappropriate")
    end

    vec(n1,n2,n3) = n1*v[1] + n2*v[2] + n3*v[3]

    l = []
    for n1 in sumRange
        for n2 in sumRange
            for n3 in zSumRange
                if n1==n2==n3==0
                    continue
                else
                    V = vec(n1,n2,n3)
                    push!(l, V)
                end
            end
        end
    end
    return l
end

function qSpaceSum(q::Vector, Δ::Vector)
    """
    Reciprocal lattice sum in d-dimensional Ewald summation

    Parameters
    ----------
    q : array_like wavevector
    Δ : array_like vector pointing between atoms in the unit cell
    Returns
    -------
    Cfar_ij : ndarray 2D array containing the reciprocal lattice sum
    """
    
end

function realSpaceSum(q::Vector, Δ::Vector)
    """
    Direct lattice sum in d-dimensional Ewald summation

    Parameters
    ----------
    q : array_like wavevector
    Δ : array_like Vector pointing between atom locations within unit cell.
    Returns
    -------
    Cfar_ij : 2D array containing the direct lattice sum
    """
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
    C_far = qSpaceSum(q, Δ)
    C_near = realSpaceSum(q, Δ)
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
