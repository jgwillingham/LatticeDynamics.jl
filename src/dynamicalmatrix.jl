


using LinearAlgebra: norm, dot




function blockMatrix(blocks::Array{Array})
    blockStack = vcat(vec(blocks)...)
    stackSize = size(blockStack)[1]
    blockSize = size(blockStack)[2]
    dimension = blockSize*size(blocks)[1]
    parts = [blockStack[i:i+dimension-1, :] for i in 1:dimension:stackSize]
    blockmatrix = hcat(parts...)
    return blockmatrix
end


# Blocks


# 𝕊_block method for shortRange interactions in bulk models
function 𝕊_block(i::Int, j::Int, q::Vector, crystal::Union{Crystal, Slab}, couplings::Array)

        A, B = couplings[i][j]
        atomᵢ = crystal.unitCell[i][1]
        atomⱼ = crystal.unitCell[j][1]
        neighborList = crystal.neighbors[atomᵢ]

        𝕊ᵢⱼ = zeros(3, 3)
        for neighbor in neighborList
                if neighbor[1] == atomⱼ
                        bondᵢⱼ = neighbor[2][1]
                        Rℓ = neighbor[2][2]
                        FCM = φ(bondᵢⱼ, A, B)
                        𝕊ᵢⱼ += FCM*exp(im*dot(q, Rℓ))
                end
        end
        return 𝕊ᵢⱼ
end

# ℂ_block method for Coulomb interactions in bulk crystals
function ℂ_block(i::Int, j::Int, q::Vector, crystal::Union{Crystal, Slab}, charges::Array)
        latticeVectors = crystal.latticeVectors
        rᵢ = crystal.unitCell[i][2]
        rⱼ = crystal.unitCell[j][2]
        Δ = dott(rⱼ - rᵢ, latticeVectors)
        Cᵢⱼ = ewald(q, Δ, crystal, charges)
end



# Self Terms


# 𝕊_selfblock method for short-range interactions in the bulk
function 𝕊_selfblock(i::Int, crystal::Union{Crystal, Slab}, couplings::Array)

        selfTerm = zeros(3,3)
        Γ = zeros(3)
        for j in eachindex(crystal.unitCell)
                selfTerm -= 𝕊_block(i, j, Γ, crystal, couplings)
        end
        return selfTerm
end




# Full dynamical matrix



# Construct the full contribution to the dynamical matrix for given interaction type
function 𝕊(q::Vector, crystal::Union{Crystal, Slab}, couplings::Array)

        atomsPerUnitCell = length(crystal.unitCell)
        blocks = Matrix{Array}(undef, (atomsPerUnitCell, atomsPerUnitCell) )
        for i in 1:atomsPerUnitCell
                for j in 1:i
                        blocks[i,j] = 𝕊_block(i, j, q, crystal, couplings)
                        if i==j
                                blocks[i,i] += 𝕊_selfblock(i, crystal, couplings)
                        else
                                blocks[j,i] = adjoint(blocks[i,j])
                        end
                end
        end
        𝕊matrix = Hermitian(blockMatrix(blocks))
        return 𝕊matrix
end


function 𝔻(q::Vector, crystal::Union{Crystal, Slab}, couplings::Array)
        𝕊ₖ = 𝕊(q, crystal, couplings)
        return 𝔻ₖ
end
