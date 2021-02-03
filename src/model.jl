


using LinearAlgebra: norm, eigvals
using Plots


function buildLine(startPoint::Vector, endPoint::Vector, pointDensity::Real)
        numPoints = pointDensity*norm(endPoint - startPoint)
        numPoints = round(Int, numPoints)
        samples = range(0, stop=1, length=numPoints)
        line = [(1-t)*startPoint + t*endPoint for t in samples]
        return line
end


function buildPath(qMarkers::Array, pointDensity::Real)
        qPathParts = Float64[0.0]
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

function getSlabCouplingArray(slab::Slab, bulkCouplingArray::Array)
        numAtoms = length(slab.unitCell)
        shift = slab.numAtomsMovedToBottom
        numBulkAtoms = length(bulkCouplingArray)
        slabCouplingArray = []
        for i in 1:numAtoms
                push!(slabCouplingArray, [])
                for j in 1:numAtoms
                        push!(slabCouplingArray[i], bulkCouplingArray[mod(i-1+shift, numBulkAtoms)+1][mod(j-1+shift, numBulkAtoms)+1])
                end
        end
        return slabCouplingArray
end


function getDispersion(qPath::Array, crystal::Union{Crystal, Slab}, couplings::Array)
        𝔻List = map(q -> 𝔻(q, crystal, couplings), qPath)
        ω²Values = map(x -> round.(x, digits=10), map(eigvals, 𝔻List))
        fValues = map( x -> .√x./(2π), ω²Values)
        return fValues
end


function plotDispersion(dispersion::Array, qPathParts::Array=[], qLabels::Array=[])
    bands = []

    for i in 1:length(dispersion[1])
        t = [set[i] for set in dispersion]
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
