


using LinearAlgebra: norm, dot


abstract type Interaction end
struct ShortRange <: Interaction end
struct Coulomb <: Interaction end




function blockMatrix(blocks::Array{Array})
    blockStack = vcat(vec(blocks)...)
    stackSize = size(blockStack)[1]
    blockSize = size(blockStack)[2]
    dimension = blockSize*size(blocks)[1]
    parts = [blockStack[i:i+dimension-1, :] for i in 1:dimension:stackSize]
    blockmatrix = hcat(parts...)
    return blockmatrix
end


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



function 𝔻_block(i::Int, j::Int, k::Vector, crystal::Crystal, couplings::Array,
                  interactionKey::ShortRange)

        A, B = couplings[i][j]
        atomᵢ = crystal.unitCell[i][1]
        atomⱼ = crystal.unitCell[j][1]
        neighborList = crystal.neighbors[atomᵢ]

        block = zeros(3, 3)
        for neighbor in neighborList
                if neighbor[1] == atomⱼ
                        bondᵢⱼ = neighbor[2]
                        FCM = φ(bondᵢⱼ, A, B)
                        block += FCM*exp(im*dot(k, bondᵢⱼ))
                end
        end
        return block
end


function 𝔻_selfblock(i::Int, crystal::Crystal, couplings::Array,
        interactionKey::ShortRange)

        selfTerm = zeros(3,3)
        Γ = zeros(3)
        for j in eachindex(crystal.unitCell)
                selfTerm -= 𝔻_block(i, j, Γ, crystal, couplings, interactionKey)
        end
        return selfTerm
end



function 𝔻_contribution(k::Vector, crystal::Crystal, couplings::Array, interactionKey::Interaction)

        atomsPerUnitCell = length(crystal.unitCell)
        blocks = Matrix{Array}(undef, (atomsPerUnitCell, atomsPerUnitCell) )
        for i in 1:atomsPerUnitCell
                for j in 1:atomsPerUnitCell
                        blocks[i,j] = 𝔻_block(i, j, k, crystal, couplings, interactionKey)
                        if i==j
                                blocks[i,j] += 𝔻_selfblock(i, crystal, couplings, interactionKey)
                        end
                end
        end
        matrix = blockMatrix(blocks)
        return matrix
end


function 𝔻(k::Vector, crystal::Crystal, couplings::Array)
        interactionKey = ShortRange()
        𝔻ₖ = 𝔻_contribution(k, crystal, couplings, interactionKey)
        return 𝔻ₖ
end
