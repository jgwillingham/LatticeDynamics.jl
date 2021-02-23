


using LinearAlgebra: norm, dot




function blockMatrix(blocks::Matrix{Matrix})
    blockStack = reduce(vcat, blocks)
    stackSize = size(blockStack)[1]
    blockSize = size(blockStack)[2]
    dimension = blockSize*size(blocks)[1]
    parts = [blockStack[i:i+dimension-1, :] for i in 1:dimension:stackSize]
    blockmatrix = reduce(hcat, parts)
    return blockmatrix
end


# Blocks


# 𝕊_block method for shortRange interactions in bulk models
function 𝕊_block(i::Int, j::Int, q::Vector{Float64}, crystal::Union{Crystal, Slab}, couplings::Array)
        e = 15.1891
        scale = e^2/(2*crystal.cellVol)
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
function ℂ_self(i::Int, crystal::Union{Crystal, Slab}, charges::Array, GList::Array, RList::Array, η::Float64)
        latticeVectors = crystal.latticeVectors
        selfTerm = zeros(3,3)
        Γ = zeros(3)
        rᵢ = crystal.cartesianUnitCell[i][2]
        for j in eachindex(crystal.unitCell)
                Zfactor = charges[j]/charges[i]
                rⱼ = crystal.cartesianUnitCell[j][2]
                Δ = rⱼ - rᵢ
                ℂᵢⱼ = ewald(Γ, Δ, crystal, GList, RList, η)
                selfTerm -= Zfactor * ℂᵢⱼ
        end
        return selfTerm
end




# Full dynamical matrix



# Construct the full contribution to the dynamical matrix from short range forces
function 𝕊(q::Vector{Float64}, crystal::Union{Crystal, Slab}, couplings::Array, atomDepth::Int)
        #atomsPerUnitCell = length(crystal.unitCell)
        blocks = Matrix{Matrix}(undef, (atomDepth, atomDepth) )
        for i in 1:atomDepth
                for j in 1:i
                        blocks[i,j] = 𝕊_block(i, j, q, crystal, couplings)
                        if i==j
                                blocks[i,j] += 𝕊_self(i, crystal, couplings)
                        else
                                blocks[j,i] = adjoint(blocks[i,j])
                        end
                end
        end
        𝕊matrix = Hermitian(blockMatrix(blocks))
        return 𝕊matrix
end



# Constructs the full coulomb contribution to the dynamical matrix
function ℂ(q::Vector, crystal::Union{Crystal, Slab}, charges::Array, sumDepth::Int, η::Float64, atomDepth::Int)
        latticeVectors = crystal.latticeVectors

        RList = getLatticeSummands(latticeVectors, sumDepth)
        GList = getLatticeSummands(crystal.reciprocalVectors, sumDepth)

        atomsPerUnitCell = length(crystal.unitCell)
        blocks = Matrix{Matrix}(undef, (atomDepth, atomDepth) )
        for i in 1:atomDepth
                rᵢ = crystal.cartesianUnitCell[i][2]
                for j in 1:i
                        rⱼ = crystal.cartesianUnitCell[j][2]
                        Δ = rⱼ - rᵢ
                        ℂᵢⱼ = ewald(q, Δ, crystal, GList, RList, η)
                        blocks[i,j] = ℂᵢⱼ
                        if i==j
                                blocks[i,i] += ℂ_self(i, crystal, charges, GList, RList, η)
                        else
                                blocks[j,i] = adjoint(blocks[i,j])
                        end
                end
        end
        ℂmatrix = Hermitian(blockMatrix(blocks))
        return ℂmatrix
end


function 𝔻(q::Vector{Float64}, crystal::Union{Crystal, Slab}, couplings::Array, atomDepth::Int=0)
        if atomDepth==0 || typeof(crystal) == Crystal{AbstractArray}
                atomDepth=length(crystal.unitCell) #the full atomDepth
        end
        𝕊ₖ = 𝕊(q, crystal, couplings, atomDepth)

        𝕄 = crystal.𝕄[1:3*atomDepth, 1:3*atomDepth]
        𝔻ₖ = Hermitian(𝕄*(𝕊ₖ)*𝕄)
        return 𝔻ₖ
end



function 𝔻(q::Vector{Float64}, crystal::Union{Crystal, Slab}, couplings::Array, charges::Array, sumDepth::Int, η::Float64, atomDepth::Int=0)
        if atomDepth==0 || typeof(crystal) == Crystal{AbstractArray}
                atomDepth=length(crystal.unitCell) #the full atomDepth
        end
        𝕊ₖ = 𝕊(q, crystal, couplings, atomDepth)
        ℂₖ = ℂ(q, crystal, charges, sumDepth, η, atomDepth)

        𝕄 = crystal.𝕄[1:3*atomDepth, 1:3*atomDepth]
        ℤ = getChargeMatrix(charges)[1:3*atomDepth, 1:3*atomDepth]
        𝔻ₖ = Hermitian( 𝕄*(𝕊ₖ + ℤ*ℂₖ*ℤ)*𝕄 )
        return 𝔻ₖ
end
