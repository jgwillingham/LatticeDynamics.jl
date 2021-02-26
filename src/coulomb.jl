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
    elseif length(latticeVectors) == 2
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

function qSpaceSum(q::Vector, Δ::Vector, crystal::Crystal, η::Float64, GList::Array)
    d = 3
    Cfar_ij = zeros(ComplexF64,3,3)
    QGList = [q+G for G in GList]

    if norm(q) > eps()
        push!(QGList, q)
    end

    for G in QGList
        qGnorm = norm(G)
        term = outer(G,G) / (qGnorm^(d-1))
        term = term * exp(-1.0im * dot(G,Δ)) #Check for using just im
        α = (d-1)/2
        x = qGnorm/(2*η)
        term = term * gamma_inc(α,x^2, 0)[2] * gamma(α) #Why the blue line?
        Cfar_ij += term
    end

    Cfar_ij = Cfar_ij * (2*sqrt(pi))^(d-1) / crystal.cellVol
end


function realSpaceSum(q::Vector, Δ::Vector, crystal::Crystal, η::Float64, RList::Array)
    Cnear_ij = zeros(ComplexF64,3,3)
    ΔRlist = [R+ Δ for R in RList]

    if norm(Δ) > sqrt(eps())
        push!(ΔRlist, Δ)
    else
        Cnear_ij += Matrix(I,3,3) * 4 / (3*sqrt(pi))
    end

    for dR in ΔRlist
        dRnorm = norm(dR)
        y = η * dRnorm
        t₁ = outer(dR,dR) / dRnorm^5
        t₁ *= (3*erfc(y) + 1/sqrt(pi) *  (6*y + 4*y^3)*exp(-y^2))
        t₂ = Matrix(I,3,3) / dRnorm^3
        t₂ *=  (erfc(y) + 2*y * exp(-y^2) / sqrt(pi))
        term = t₁ - t₂
        term *= exp(1.0im * dot(q,dR-Δ) )
        Cnear_ij += term
    end
    return -1*Cnear_ij
end


function ewald(q::Vector, Δ::Vector, crystal::Crystal, GList::Array, RList::Array, η::Float64)
    C_far = qSpaceSum(q, Δ, crystal, η, GList)
    C_near = realSpaceSum(q, Δ, crystal, η, RList)
    C_ij = C_far + C_near
end


function ewald(q::Vector, Δ::Vector, crystal::Slab, GList::Array, RList::Array, η::Float64)
    Δₚ, Δₙ = projectVector(Δ, crystal.surfaceNormal)

    if norm(Δₙ) > sqrt(eps())
        C_ij = differentPlaneSum(q, Δₚ, Δₙ, crystal, GList)
    else
        C_ij = samePlaneSumDeWette(q, Δ, crystal, GList, RList)
    end
    C_ij
end


function differentPlaneSum(q, Δparallel::Vector, Δnormal::Vector, crystal::Slab, GList::Array{Vector})
    C_ij = zeros(ComplexF64,3,3)
    qGList = [q+G for G in GList] #Check where to find

    if norm(q) > sqrt(eps())
        push!(qGList, q)
    end
    n = crystal.surfaceNormal
    sign = sign(dot(n, Δnormal))

    for qG in qGList
        qGnorm = norm(qG)
        Δₙ = norm(Δnormal)

        t₁ = outer(qG, qG)/qGnorm
        t₂ = 1.0im * outer(qG,n)
        t₂ += 1.0im*outer(n,qG)
        t₃ = qGnorm*outer(n,n)
        t = (t₁ - sign*t₂ - t₃)*exp(-1.0im*dot(qG,Δparallel))
        t = t * exp( - Δₙ * qGnorm)

        C_ij = C_ij + t
    end
    C_ij = C_ij * (2*pi / crystal.meshArea)
end


# Equivalent of _DeWette in coulomb.py
function samePlaneSumDeWette(q::Vector, Δ::Vector, crystal::Slab, GList::Array{Vector}, RList::Array{Vector})
    id_xy0 = [1 0 0;
              0 1 0;
              0 0 0]
    Ec = [1 0 0;
          0 1 0;
          0 0 -2]
    a = norm(crystal.meshPrimitives[1])
    Δ = Δ/a
    q = q*a

    C_ij = 2*pi/(crystal.meshArea/a) * Ec
    if norm(Δ) < sqrt(eps())
        C_ij -= 2*pi/3 * Ec
    end
    Cfar_ij = zeros(3,3)
    qGList = [q+a*G for G in GList]

    for qG in qGList
        qGnorm = norm(qG)
        arg = qGnorm^2 / (4*pi)
        ϕ = exp(-1.0im * dot(qG,Δ))
        t₁ = gamma_inc(1/2,arg,0)[1] * (2* outer(qG,qG)/qGnorm^2 - id_xy0)
        t₂ = Matrix(I,3,3)/(-2) * gamma_inc(-0.5,arg,0)[1]
        t = qGnorm*ϕ*(t₁+t₂)
        Cfar_ij = Cfar_ij + t
    end
    Cfar_ij *= -sqrt(pi)/(crystal.meshArea/a)

    Cnear_ij = zeros(ComplexF64,3,3)

    DeltaRList = [(R/a + Δ) for R in RList]
    if norm(Δ) > sqrt(eps())
        push!(DeltaRList, Δ) # include R=0 term when non-singular
    end

    for dR in DeltaRList
        norm = norm(dR)
        arg = pi*norm^2
        ϕ = exp(1im* dot(q,(dR-Δ) )) # only lattice vector appears in phase (questionable: Lucas)

        t₁ = gammainc(5/2, arg) * ( 2*outer(dR, dR)/norm^2 - id_xy0 )
        t₂ = Matrix(I,3,3)/2 * gamma_inc(3/2, arg,0)[1]
        t = (ϕ/norm^3) * (t₁ + t₂)

        Cnear_ij = Cnear_ij + t
    end
    # got rid of minus sign below
    Cnear_ij *= 2/sqrt(pi)

    C_ij = C_ij + (Cnear_ij + Cfar_ij)*Ec

    C_ij = C_ij / (a^3)
    return -1*C_ij
end
