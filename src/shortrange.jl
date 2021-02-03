


using LinearAlgebra: norm, transpose

# outer product of two vectors
function outer(v1::Vector, v2::Vector)
        return v1.*transpose(v2)
end


# short range radial force constant matrices
function φ(bondᵢⱼ::Vector, A::Real, B::Real)
        bondLength = norm(bondᵢⱼ)
        ϕ = (A-B)*outer(bondᵢⱼ, bondᵢⱼ) / bondLength^2
        ϕ += B * Matrix(I,3,3)
        return ϕ
end
