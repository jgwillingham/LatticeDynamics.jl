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

function qSpaceSum(q::Vector, Δ::Vector,η::Float64, d::Int32, GList::Array{Vector})
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
    Cfar_ij = zeros(3,3){ComplexF64}
    QGList = [q+G for G in GList] 

    if norm(q) > eps()
        push!(QGList, q)
    end

    for G in QGList
        norm = norm(G)
        term = outer(G,G) / (norm^(d-1)) 
        term = term * exp(-1.0im * dot(G,Δ)) #Check for using just im
        α = (d-1)/2
        x = norm/(2*η)
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
        t₂ *=  (erfc(y) + 2*y * exp(-y^2) / sqrt(pi))
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
    """
    """
    Δₚ, Δₙ = SELF.LATTICE.PROJECTVECTOR(Δ) #Check where to find

    if norm(Δₙ) > sqrt(eps())
        C_ij = differentPlaneSum(q,Δₚ,Δₙ)
    else
        C_ij = samePlaneSumDeWette(q,Δ,crystal)
    end
    C_ij
end


function differentPlaneSum(q, Δparallel::Vector, Δnormal::Vector, crystal::Slab)
    """
    """
    C_ij = zeros(3,3){ComplexF64}
    qGList = [q+G for G in SELF.GLIST] #Check where to find

    if norm(q) > sqrt(eps())
        push!(qGList, q)
    end
    n = SELF.LATTICE.SURFACENORMAL #Check
    sign = sign(n * Δnormal) 

    for qG in qGList
        qGnorm = norm(qG)
        Δₙ = norm(Δnormal) #Check

        t₁ = outer(qG, qG)/qGnorm
        t₂ = 1im * outer(qG,n)
        t₂ += 1im*outer(n,qG)
        t₃ = qGnorm*outer(n,n)
        t = (t₁ - sign*t₂ - t₃)*exp(-1im*dot(qG,Δparallel)) 
        t = t * exp( - Δₙ * qGnorm)
        
        C_ij = C_ij + t
    end
    C_ij = C_ij * (2*pi / SELF.CELLVOL) #Check
end


# Equivalent of _DeWette in coulomb.py
function samePlaneSumDeWette(q::Vector, Δ::Vector, crystal::Slab)
    """
    """
    id_xy0 = [[1 0 0;
               0 1 0;
               0 0 0]]
    Ec = [[1 0 0;
           0 1 0;
           0 0 -2]]
    a = norm(SELF.LATTICE.MESHPRIMITIVES[1]) #Check where to find
    Δ = Δ/a
    q = q*a

    C_ij = 2*pi/(SELF.CELLVOL/a) * Ec 
    if norm(Δ) < sqrt(eps()) #If there are problems with the value being too high, use eps()
        C_ij -= 2*pi/3 * Ec
    end
    Cfar_ij = zeros(3,3)
    qGList = [q+a*G for G in SELF.GLIST] #Check where to find

    for qG in qGList
        qGnorm = norm(qG)
        arg = qGnorm^2 / (4*pi)
        ϕ = exp(-1im * dot(qG,Δ))
        t₁ = gamma_inc(1/2,arg)[1] * (2* outer(qG,qG)/qGnorm^2 - id_xy0) #Check if float casting necessary on gamma_inc
        t₂ = Matrix(I,3,3)/(-2) * gamma_inc(-0.5,arg)[1]
        t = qGnorm*ϕ*(t₁+t₂)
        Cfar_ij = Cfar_ij + t
    end
    Cfar_ij *= -sqrt(pi)/(SELF.CELLVOL/a) #Check where to find

    Cnear_ij = zeros(3,3){ComplexF64}
        
    DeltaRList = [(R/a + Δ) for R in SELF.RLIST] #Check where to find
    if norm(Δ) > sqrt(eps())
        push!(DeltaRList, Δ) # include R=0 term when non-singular
    end

    for dR in DeltaRList
        norm = norm(dR)
        arg = pi*norm^2
        ϕ = exp(1im* dot(q,(dR-Δ) )) # only lattice vector appears in phase (questionable: Lucas)
        
        t₁ = gammainc(5/2, arg) * ( 2*outer(dR, dR)/norm^2 - id_xy0 )  
        t₂ = Matrix(I,3,3)/2 * gamma_inc(3/2, arg)[1]
        t = (ϕ/norm^3) * (t₁ + t₂) # times Ec?

        Cnear_ij = Cnear_ij + t
    end
    # got rid of minus sign below
    Cnear_ij *= 2/sqrt(pi)

    C_ij = C_ij + (Cnear_ij + Cfar_ij)*Ec  #Originally was @ Ec, I think this works

    C_ij = C_ij / (a^3) # scaling
    return -1*C_ij
end
