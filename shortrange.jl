module Model

using Plots
using LinearAlgebra

export Crystal, getNeighbors!, buildPath, getDispersion, plotDispersion


function dott(x::Array, y::Array)
        res = sum(x[i].*y[i] for i in range(1, stop=length(x)))
        return res
end



function blockMatrix(blocks::Array{Array})
    blockStack = vcat(vec(blocks)...)
    stackSize = size(blockStack)[1]
    blockSize = size(blockStack)[2]
    dimension = blockSize*size(blocks)[1]
    parts = [blockStack[i:i+dimension-1, :] for i in 1:dimension:stackSize]
    blockmatrix = hcat(parts...)
    return blockmatrix
end


mutable struct Crystal
        unitCell::Array
        latticeVectors::Array
        neighbors::Dict
        function Crystal(unitCell, latticeVectors)
                cartesian_unitCell = [[atom[1], dott(atom[2],latticeVectors)]
                                        for atom in unitCell]
                new(cartesian_unitCell, latticeVectors, Dict())
        end
end


function getNeighbors!(lattice::Crystal, threshold::Float64, searchWidth::Int=2)
        searchRange = range(-searchWidth, stop=searchWidth)
        neighbors = Dict()
        for atomᵢ in lattice.unitCell
                neighbors[atomᵢ[1]] = []
                rᵢ = atomᵢ[2]
                for atomⱼ in lattice.unitCell
                        xⱼ = atomⱼ[2]
                        for n₁ in searchRange
                                for n₂ in searchRange
                                        for n₃ in searchRange
                                                Rℓ = dott([n₁, n₂, n₃], lattice.latticeVectors)
                                                rⱼ  = Rℓ + xⱼ
                                                distance = norm(rⱼ - rᵢ)
                                                if distance < threshold && distance !=0.0
                                                        append!(neighbors[atomᵢ[1]], [[atomⱼ[1], rⱼ-rᵢ]])
                                                end
                                        end
                                end
                        end
                end
        end
        lattice.neighbors = neighbors
end



function φ(bondᵢⱼ::Vector, A::Float64, B::Float64)
        ϕ = zeros(3, 3)
        bondLength = norm(bondᵢⱼ)
        for xμ in range(1, stop=3)
                for xν in range(1, stop=3)
                        ϕ[xμ, xν] = (A-B)*bondᵢⱼ[xμ]*bondᵢⱼ[xν] / bondLength^2
                        if xμ==xν
                                ϕ[xμ, xν] += B
                        end
                end
        end
        return ϕ
end



function 𝔻ₖ_block(k::Vector, i::Int, j::Int,
                  crystal::Crystal, couplings::Array)

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


function 𝔻ₖ_selfblock(i::Int, crystal::Crystal, couplings::Array)
        selfTerm = zeros(3,3)
        Γ = zeros(3)
        for j in range(1, stop=length(crystal.unitCell))
                selfTerm -= 𝔻ₖ_block(Γ, i, j, crystal, couplings)
        end
        return selfTerm
end



function 𝔻ₖ(k::Vector, crystal::Crystal, couplings::Array)
        atomsPerUnitCell = length(crystal.unitCell)
        blocks = Matrix{Array}(undef, (atomsPerUnitCell, atomsPerUnitCell) )
        for i in 1:atomsPerUnitCell
                for j in 1:atomsPerUnitCell
                        blocks[i,j] = 𝔻ₖ_block(k, i, j, crystal, couplings)
                        if i==j
                                blocks[i,j] += 𝔻ₖ_selfblock(i, crystal, couplings)
                        end
                end
        end
        matrix = blockMatrix(blocks)
        return matrix
end


function getDispersion(qPath::Array, crystal::Crystal, couplings::Array)
        eigenValues = []
        for q in qPath
                𝔻 = 𝔻ₖ(q, crystal, couplings)
                vals = map( x -> round(x, digits=10) , eigvals(𝔻))
                freq = .√vals./(2π)
                append!(eigenValues, [freq])
        end
        return eigenValues
end


function buildLine(startPoint::Vector, endPoint::Vector, pointDensity::Real)
        numPoints = pointDensity*norm(endPoint - startPoint)
        numPoints = round(Int, numPoints)
        samples = range(0, stop=1, length=numPoints)
        line = [(1-t)*startPoint + t*endPoint for t in samples]
        return line
end


function buildPath(qMarkers::Array, pointDensity::Real)
        qPathParts = [0.0]
        firstLine = buildLine(qMarkers[1], qMarkers[2], pointDensity)
        qPath = firstLine
        append!(qPathParts, length(firstLine))
        for i in 2:length(qMarkers)-1
                startPoint = qMarkers[i]
                endPoint = qMarkers[i+1]
                line = buildLine(startPoint, endPoint, pointDensity)[2:end]
                append!(qPath, line)
                append!(qPathParts, length(qPath))
        end
        return qPath, qPathParts
end


function plotDispersion(disp::Array, qPathParts::Array=[], qLabels::Array=[])
    bands = []

    for i in 1:length(disp[1])
        t = [set[i] for set in disp]
        t = map(real, t)
        append!(bands, [t])
    end
    plot(bands, linewidth=2, xticks=(qPathParts, qLabels), legend=false, xtickfont=(13), size=(800, 400))
    plot!(qPathParts, seriestype=:vline, color=:black, linealpha=0.35)
    xlims!((0.0, qPathParts[end]))
    ylims!(0, Inf)
    title!("Phonon Dispersion")
    ylabel!("Frequency (THz)")
end

end # module
