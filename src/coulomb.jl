

# this will contain all functions for implementing Ewald summation

function getChargeMatrix(charges::Array)
    "Get charge for matrix Z, returning a diagonal
    matrix containing the charges"
    e = 15.1891
    Z = Diagonal(repeat(charges, inner=[3])) .* e
end


# Returns a list of lattice (real or reciprocal) vectors which are to be summed over.
# Equivalent of _buildList in coulomb.py
function getLatticeSummands(crystal::Crystal, sumDepth::Int)
    "Builds and returns a list of vectors to be summed over in Ewald summation"
    sumRange = -sumDepth:1:sumDepth
    a₁, a₂, a₃ = crystal.latticeVectors
    b₁, b₂, b₃ = crystal.reciprocalVectors

    latVec(n1,n2,n3) = n1*a₁ + n2*a₂ + n3*a₃
    recipVec(n1,n2,n3) = n1*b₁ + n2*b₂ + n3*b₃

    realVecs = []
    recipVecs = []
    for n1 in sumRange
        for n2 in sumRange
            for n3 in sumRange
                if n1==n2==n3==0
                    continue
                else
                    R = latVec(n1,n2,n3)
                    push!(realVecs, R)
                    G = recipVec(n1,n2,n3)
                    push!(recipVecs, G)
                end
            end
        end
    end
    return realVecs, recipVecs
end


function getLatticeSummands(slab::Slab, sumDepth::Int)
    "Builds and returns a list of vectors to be summed over in Ewald summation"
    sumRange = -sumDepth:1:sumDepth
    a₁, a₂ = slab.meshPrimitives
    b₁, b₂ = slab.meshReciprocals

    latVec(n1,n2) = n1*a₁ + n2*a₂
    recipVec(n1,n2) = n1*b₁ + n2*b₂

    realVecs = []
    recipVecs = []
    for n1 in sumRange
        for n2 in sumRange
            if n1==n2==0
                continue
            else
                R = latVec(n1,n2)
                push!(realVecs, R)
                G = recipVec(n1,n2)
                push!(recipVecs, G)
            end
        end
    end
    return realVecs, recipVecs
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
        term = term * sf_gamma_inc(α,x^2)
        Cfar_ij += term
    end

    Cfar_ij = Cfar_ij * (2*sqrt(pi))^(d-1) / crystal.cellVol
end


function qSpaceSum(q::Vector, Δ::Vector, slab::Slab, η::Float64, GList::Array)
    d = 2
    Cfar_ij = zeros(ComplexF64,3,3)
    QGList = [q + G for G in GList]

    if norm(q) > eps()
        push!(QGList, q)
    end

    for G in QGList
        qGnorm = norm(G)
        term = outer(G,G) / (qGnorm^(d-1))
        term = term * exp(-1.0im * dot(G,Δ)) #Check for using just im
        α = (d-1)/2
        x = qGnorm/(2*η)
        term = term * sf_gamma_inc(α,x^2)
        Cfar_ij += term
    end

    Cfar_ij = Cfar_ij * (2*sqrt(pi))^(d-1) / slab.meshArea
end


function realSpaceSum(q::Vector, Δ::Vector, η::Float64, RList::Array)
    Cnear_ij = zeros(ComplexF64,3,3)
    ΔRlist = [R + Δ for R in RList]

    if norm(Δ) > sqrt(eps())
        push!(ΔRlist, Δ)
    #else
        # Cnear_ij += Matrix(I,3,3) * 4 / (3*sqrt(pi))
    end

    for dR in ΔRlist
        dRnorm = norm(dR)
        y = η * dRnorm
        t₁ = outer(dR,dR) / dRnorm^5
        t₁ *= (3*sf_erfc(y) + 1/sqrt(pi) *  (6*y + 4*y^3)*exp(-y^2))
        t₂ = Matrix(I,3,3) / dRnorm^3
        t₂ *=  (sf_erfc(y) + 2*y * exp(-y^2) / sqrt(pi))
        term = t₁ - t₂
        term *= exp(1.0im * dot(q,dR-Δ) )
        Cnear_ij += term
    end
    return -1*Cnear_ij
end


function ewald(q::Vector, Δ::Vector, crystal::Crystal, GList::Array, RList::Array, η::Float64)
    C_far = qSpaceSum(q, Δ, crystal, η, GList)
    C_near = realSpaceSum(q, Δ, η, RList)
    C_ij = C_far + C_near
end


function ewald(q::Vector, Δ::Vector, crystal::Slab, GList::Array, RList::Array, η::Float64)
    Δₚ, Δₙ = projectVector(Δ, crystal.surfaceNormal)

    if norm(Δₙ) > sqrt(eps())
        C_ij = differentPlaneSum(q, Δₚ, Δₙ, crystal, GList)
    else
        #C_ij = samePlaneSumDeWette(q, Δ, crystal, GList, RList)
        C_far = qSpaceSum(q, Δ, crystal, η, GList)
        C_near = realSpaceSum(q, Δ, η, RList)
        C_ij = C_far + C_near
    end
    C_ij
end


function differentPlaneSum(q::Vector, Δparallel::Vector, Δnormal::Vector, crystal::Slab, GList::Array)
    C_ij = zeros(ComplexF64,3,3)
    qGList = [q+G for G in GList]

    if norm(q) > sqrt(eps())
        push!(qGList, q)
    end
    n = crystal.surfaceNormal
    sgn = sign(dot(n, Δnormal))

    for qG in qGList
        qGnorm = norm(qG)
        Δₙ = norm(Δnormal)

        t₁ = outer(qG, qG)/qGnorm
        t₂ = 1.0im * outer(qG,n)
        t₂ += 1.0im*outer(n,qG)
        t₃ = qGnorm*outer(n,n)
        t = (t₁ - sgn*t₂ - t₃)*exp(-1.0im*dot(qG,Δparallel))
        t = t * exp( - Δₙ * qGnorm)

        C_ij = C_ij + t
    end
    C_ij = C_ij * (2*pi / crystal.meshArea)
end


# Equivalent of _DeWette in coulomb.py
function samePlaneSumDeWette(q::Vector, Δ::Vector, crystal::Slab, GList::Array, RList::Array)
    id_xy0 = [1. 0. 0.;
              0. 1. 0.;
              0. 0. 0.]
    Ec = [1. 0. 0.;
          0. 1. 0.;
          0. 0. -2.]
    a = norm(crystal.meshPrimitives[1])
    Δ = Δ/a
    q = q*a

    C_ij = 2*pi/(crystal.meshArea/a) * Ec
    if norm(Δ) < sqrt(eps())
        C_ij -= 2*pi/3 * Ec
    end
    Cfar_ij = zeros(ComplexF64,3,3)
    qGList = [q+a*G for G in GList]

    for qG in qGList
        qGnorm = norm(qG)
        arg = qGnorm^2 / (4*pi)
        ϕ = exp(-1.0im * dot(qG,Δ))
        t₁ = sf_gamma_inc(1/2, arg) * (2* outer(qG,qG)/qGnorm^2 - id_xy0)
        t₂ = Matrix(I,3,3)/(-2) * sf_gamma_inc(-1/2, arg)
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
        dRnorm = norm(dR)
        arg = pi*dRnorm^2
        ϕ = exp(1.0im* dot(q,(dR-Δ) )) # only lattice vector appears in phase (questionable: Lucas)

        t₁ = sf_gamma_inc(5/2, arg) * ( 2*outer(dR, dR)/dRnorm^2 - id_xy0 )
        t₂ = Matrix(I,3,3)/2 * sf_gamma_inc(3/2, arg)
        t = (ϕ/dRnorm^3) * (t₁ + t₂)

        Cnear_ij = Cnear_ij + t
    end
    # got rid of minus sign below
    Cnear_ij *= 2/sqrt(pi)

    C_ij = C_ij + (Cnear_ij + Cfar_ij)*Ec

    C_ij = C_ij / (a^3)
    return -1*C_ij
end
