


function getChargeMatrix(charges::Array)
    "Get charge for matrix Z, returning a diagonal
    matrix containing the charges"
    e = 15.1891
    Z = Diagonal(repeat(charges, inner=[3])) .* e
end


@inline function RSumTerm(ΔR::Vector, η::Float64)
    ΔRnorm = norm(ΔR)
    y = η * ΔRnorm
    C = (outer(ΔR, ΔR) / ΔRnorm^5) * (3*sf_erfc(y) + 1/√pi *  (6*y + 4*y^3)*exp(-y^2))
    C -=  (Matrix(I,3,3)/ ΔRnorm^3) * (sf_erfc(y) + 2*y * exp(-y^2) / √π )
    return -1*C
end


@inline function GSumTerm(qG::Vector, η::Float64)
    qGnorm = norm(qG)
    x = qGnorm/(2*η)
    C = outer(qG,qG) / qGnorm^2 * exp(-x^2)
    return C
end


@inline function slabGSumTerm(qG::Vector, η::Float64)
    qGnorm = norm(qG)
    x = qGnorm/(2*η)
    C = outer(qG,qG) / qGnorm * sf_erfc(x)
    return C
end

#bulk ewald
function ewald(q::Vector, Δ::Vector, crystal::Crystal, sumDepth::Int, η::Float64)
    sumRange = -sumDepth:1:sumDepth
    a₁, a₂, a₃ = crystal.latticeVectors
    b₁, b₂, b₃ = crystal.reciprocalVectors

    sumList = [[n1,n2,n3] for n1 in sumRange for n2 in sumRange for n3 in sumRange]
    filter!(e->e!=zeros(3), sumList)

    C_ij = zeros(ComplexF64, 3, 3)
    for integers in sumList
        n1, n2, n3 = integers

        Rℓ = n1*a₁ + n2*a₂ + n3*a₃
        ΔR = Δ + Rℓ
        C_ij += RSumTerm(ΔR, η) * exp(im*dot(q, Rℓ) )

        G = n1*b₁ + n2*b₂ + n3*b₃
        qG = q + G
        C_ij += 4π/crystal.cellVol * GSumTerm(qG, η) * exp(-im*dot(qG, Δ))
    end
    # handle zero terms
    if norm(Δ) > √eps()
        C_ij -= RSumTerm(Δ, η)
    end
    if norm(q) > √eps()
        C_ij += 4π/crystal.cellVol * GSumTerm(q, η) * exp(-im*dot(q, Δ))
    end

    return C_ij
end


# slab Ewald
function ewald(q::Vector, Δ::Vector, slab::Slab, sumDepth::Int, η::Float64)
    sumRange = -sumDepth:1:sumDepth
    a₁, a₂ = slab.meshPrimitives
    b₁, b₂ = slab.meshReciprocals

    sumList = [[n1,n2] for n1 in sumRange for n2 in sumRange]
    filter!(e->e!=zeros(2), sumList)

    n = slab.surfaceNormal
    Δparallel, Δnormal = projectVector(Δ, n)

    C_ij = zeros(ComplexF64, 3, 3)
    if norm(Δnormal) > √eps() # atoms in different planes --> all in reciprocal space
        for integers in sumList
            n1, n2 = integers
            G = n1*b₁ + n2*b₂
            qG = q + G
            C_ij += 2π/slab.meshArea * differentPlaneSumTerm(qG, Δnormal, n) * exp(-im*dot(qG, Δparallel))
        end
        if norm(q) > √eps()
            C_ij += 2π/slab.meshArea * differentPlaneSumTerm(q, Δnormal, n) * exp(-im*dot(q, Δparallel))
        end
    else # atoms in same plane  --> Ewald method
        for integers in sumList
            n1, n2 = integers
            Rℓ = n1*a₁ + n2*a₂
            ΔR = Δ + Rℓ
            C_ij += RSumTerm(ΔR, η) * exp(im*dot(q, Rℓ) )
            G = n1*b₁ + n2*b₂
            qG = q + G
            C_ij += 2π/slab.meshArea * slabGSumTerm(qG, η) * exp(-im*dot(qG, Δ))
        end
        if norm(q) > √eps() # include G=0 term when q != 0
            C_ij += 2π/slab.meshArea * slabGSumTerm(q, η) * exp(-im*dot(q, Δ))
        end
        if norm(Δ) > √eps() # include Rℓ=0 term when Δ != 0 (i.e. not the same atom)
            C_ij += RSumTerm(Δ, η)
        else # when |Δ|=0, replace Rℓ=0 term with 4/(3√π)I
            C_ij += 4/(3*√π) * Matrix(I,3,3)
        end

        # DeWette calculation: REQUIRES A SLAB WITH SURFACE NORMAL IN Z DIRECTION

        # a = norm(a₁) # normalization factor
        # for integers in sumList
        #     n1, n2 = integers
        #     Rℓ = n1*a₁ + n2*a₂
        #     σ = (Δ + Rℓ)/a
        #     C_ij += -2/√π * RSumTermDW(σ) * exp(im*dot(q, Rℓ))
        #     G = n1*b₁ + n2*b₂
        #     h = a*(q+G)
        #     C_ij += √π/ (slab.meshArea/a) * GSumTermDW(h) * exp(-im*dot(q+G, Δ))
        # end
        # C_ij += 2π/(slab.meshArea/a) * I
        # if norm(Δ) > √eps()
        #     C_ij += -2/√π * RSumTermDW(Δ/a)
        # else
        #     C_ij += 2π/3 * I
        # end
        # C_ij *= (1/a^3) * [1. 0. 0.; 0. 1. 0.; 0. 0. -2.] # this is the Ec matrix
    end

    return C_ij
end




@inline function differentPlaneSumTerm(qG::Vector, Δnormal::Vector, surfaceNormal::Vector)
    n = surfaceNormal
    sgn = sign(dot(n, Δnormal))
    qGnorm = norm(qG)
    Δnorm = norm(Δnormal)

    C = outer(qG, qG)/qGnorm - im*sgn*(outer(qG,n) + outer(n,qG)) - qGnorm*outer(n,n)
    C *= exp(-Δnorm*qGnorm)
    return C
end


@inline function RSumTermDW(σ::Vector)
    id_xy0 = [1. 0. 0.; 0. 1. 0.; 0. 0. 0.]
    σnorm = norm(σ)
    x = π*σnorm^2
    C = sf_gamma_inc(5/2, x) * ( 2*outer(σ, σ)/σnorm^2 - id_xy0 )
    C += I/2 * sf_gamma_inc(3/2, x)
    C *= 1/σnorm^3
    return C
end


@inline function GSumTermDW(h::Vector)
    id_xy0 = [1. 0. 0.; 0. 1. 0.; 0. 0. 0.]
    hnorm = norm(h)
    x = hnorm^2/(4π)
    C  = sf_gamma_inc(1/2, x) * ( 2*outer(h, h)/hnorm^2 - id_xy0 )
    C += -I/2 * sf_gamma_inc(-1/2, x)
    C *= hnorm
    return C
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
