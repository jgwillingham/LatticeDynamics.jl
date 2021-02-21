using LinearAlgebra: norm
using SpecialFunctions

# this will contain all functions for implementing Ewald summation

function getChargeMatrix(charges::Array)
    "Get charge for matrix Z, returning a diagonal
    matrix containing the charges"
    e = 15.1891
    Z = Diagonal(repeat(charges, inner=[3])) .* e
end


# Returns a list of lattice (real or reciprocal) vectors which are to be summed over. 
# Equivalent of _buildList in coulomb.py
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
    d = SELF.DIM #check where to get dim
    Cfar_ij = zeros(3,3){ComplexF64}
    QGList = [q+G for G in SELF.GLIST] #Check where to get GList + dim q+G

    if norm(q) > sqrt(eps())
        push!(QGList, q)
    end

    for G in QGList
        norm = norm(G)
        term = outer(G,G) / (norm^(d-1)) #check if broadcasting needed
        term = term * exp(-1im * dot(G,Δ)) #check dot/inner
        α = (d-1)/2
        x = norm/(2*SELF.ETA)
        term = term * gamma_inc(α,x^2)[2] * gamma(α) #Why the blue line?
        Cfar_ij += term
    end

    Cfar_ij = Cfar_ij * (2*sqrt(pi))^(d-1) / SELF.CELLVOL
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
    Cnear_ij : 2D array containing the direct lattice sum
    """
    Cnear_ij = zeros(3,3){ComplexF64}
    ΔRlist = [R+ Δ for R in SELF.RLIST]  #Check RLIST

    if norm(Δ) > 10^(-9)
        push!(ΔRList,Δ)
    else
        Cnear_ij += Matrix(I,3,3) * 4 / (3*sqrt(pi))
    end

    for dR in ΔRlist
        norm = norm(dR)
        y = SELF.ETA * norm   #Check eta
        t₁ = outer(dR,dR) / norm^5
        t₁ *= (3*erfc(y) + 1/sqrt(pi) *  (6*y + 4*y^3)*exp(-y^2))
        t₂ = Matrix(I,3,3) / norm^3
        t₂ *=  ( erfc(y) + 2*y * exp(-y^2) / sqrt(pi) )
        term = t1 - t2
        term *= exp(1im * dot(q,dR-Δ) )
        Cnear_ij += term
    end
    return -1*Cnear_ij
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
