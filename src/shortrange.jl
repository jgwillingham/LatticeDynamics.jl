


using LinearAlgebra: norm

# short range radial force constant matrices
function φ(bondᵢⱼ::Vector, A::Real, B::Real)
        ϕ = zeros(3, 3)
        bondLength = norm(bondᵢⱼ)
        for xμ in 1:3
                for xν in 1:3
                        ϕ[xμ, xν] = (A-B)*bondᵢⱼ[xμ]*bondᵢⱼ[xν] / bondLength^2
                        if xμ == xν
                                ϕ[xμ, xν] += B
                        end
                end
        end
        return ϕ
end
