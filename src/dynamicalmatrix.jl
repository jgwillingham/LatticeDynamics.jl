


using LinearAlgebra: norm, dot
using BlockArrays: BlockArray, setblock!




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
        e = 15.1891
        a₁, a₂, a₃ = crystal.latticeVectors
        vol = abs(dot(a₁, cross(a₂, a₃)))
        scale = e^2/(2*vol)
        A, B = couplings[i][j]
        A *= scale
        B *= scale
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


# SELF TERMS


# 𝕊_self method for short-range interactions in the bulk
function 𝕊_self(i::Int, crystal::Union{Crystal, Slab}, couplings::Array)

        selfTerm = zeros(3,3)
        Γ = zeros(3)
        for j in eachindex(crystal.unitCell)
                selfTerm -= 𝕊_block(i, j, Γ, crystal, couplings)
        end
        return selfTerm
end


# ℂ_self for coulomb self interaction
function ℂ_self(i::Int, crystal::Union{Crystal, Slab}, charges::Array)
        latticeVectors = crystal.LatticeVectors
        selfTerm = zeros(3,3)
        Γ = zeros(3)
        rᵢ = dott(crystal.unitCell[i][2], latticeVectors)
        for j in eachindex(crystal.unitCell)
                Zfactor = charges[j]/charges[i]
                rⱼ = dott(crystal.unitCell[j][2], latticeVectors)
                Δ = rⱼ - rᵢ
                ℂᵢⱼ = ewald(Γ, Δ, crystal, charges)
                selfTerm -= Zfactor * ℂᵢⱼ
        end
        return selfTerm
end




# Full dynamical matrix



# Construct the full contribution to the dynamical matrix from short range forces
function 𝕊(q::Vector, crystal::Union{Crystal, Slab}, couplings::Array)
        atomsPerUnitCell = length(crystal.unitCell)
        blocks = Matrix{Array}(undef, (atomsPerUnitCell, atomsPerUnitCell) )
        #blocks = BlockArray{Complex}(undef, repeat([3], atomsPerUnitCell), repeat([3], atomsPerUnitCell) )
        for i in 1:atomsPerUnitCell
                for j in 1:i
                        blocks[i,j] = 𝕊_block(i, j, q, crystal, couplings)
                        #blockᵢⱼ = 𝕊_block(i, j, q, crystal, couplings)
                        if i==j
                                #blockᵢⱼ += 𝕊_self(i, crystal, couplings)
                                blocks[i,j] += 𝕊_self(i, crystal, couplings)
                                #setblock!(blocks, blockᵢⱼ, i, j)
                        else
                                #blockⱼᵢ = adjoint(blockᵢⱼ)
                                blocks[j,i] = adjoint(blocks[i,j])
                                #setblock!(blocks, blockᵢⱼ, i, j)
                                #setblock!(blocks, blockⱼᵢ, j, i)
                        end
                end
        end
        #𝕊matrix = Hermitian(Array(blocks))
        𝕊matrix = Hermitian(blockMatrix(blocks))
        return 𝕊matrix
end



# Constructs the full coulomb contribution to the dynamical matrix
function ℂ(q::Vector, crystal::Union{Crystal, Slab}, charges::Array)
        latticeVectors = crystal.latticeVectors
        atomsPerUnitCell = length(crystal.unitCell)
        blocks = Matrix{Array}(undef, (atomsPerUnitCell, atomsPerUnitCell) )
        for i in 1:atomsPerUnitCell
                rᵢ = dott(crystal.unitCell[i][2], latticeVectors)
                for j in 1:i
                        rⱼ = dott(crystal.unitCell[j][2], latticeVectors)
                        Δ = rⱼ - rᵢ
                        ℂᵢⱼ = ewald(q, Δ, crystal, charges)
                        blocks[i,j] = ℂᵢⱼ
                        if i==j
                                blocks[i,i] += ℂ_self(i, crystal, charges)
                        else
                                blocks[j,i] = adjoint(blocks[i,j])
                        end
                end
        end
        ℂmatrix = Hermitian(blockMatrix(blocks))
        return ℂmatrix
end


function 𝕕(q::Vector, crystal::Union{Crystal, Slab}, couplings::Array)
        𝕊ₖ = 𝕊(q, crystal, couplings)
        # ℂₖ = ℂ(q, crystal, charges)
        𝕕ₖ = 𝕊ₖ #+ ℂₖ
        return 𝕕ₖ
end
