


function buildLine(startPoint::Vector, endPoint::Vector, pointDensity::Real)
        numPoints = pointDensity*norm(endPoint - startPoint)
        numPoints = round(Int, numPoints)
        samples = range(0., stop=1., length=numPoints)
        line = [(1-t)*startPoint + t*endPoint for t in samples]
        return line
end


function buildPath(qMarkers::Array, pointDensity::Real)
        qPathParts = Float64[1.0]
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


function getSlabCharges(slab::Slab, charges::Array)
        slabCharges = repeat(charges, slab.numCells)
        rolledSlabCharges = circshift(slabCharges, slab.numAtomsMovedToBottom)
        return rolledSlabCharges
end


function getDispersion(qPath::Array, crystal::Union{Crystal, Slab}, couplings::Array)
        𝔻List = map(q -> 𝔻(q, crystal, couplings), qPath)
        eigeninfo = eigen.(𝔻List)
        ω²Values = [round.(eig.values, digits=5) for eig in eigeninfo]
        fValues = map( x -> .√Complex.(x)./(2π), ω²Values)
        meVDispersion = 4.13567 .*fValues # convert THz to meV

        normalModes = [eig.vectors for eig in eigeninfo]

        return meVDispersion, normalModes
end


function getDispersion(qPath::Array, crystal::Union{Crystal, Slab}, couplings::Array, charges::Array, sumDepth::Int=5, η=nothing, showProgress::Bool=true)
        if η == nothing
                η = 4*crystal.cellVol^(-1/3) # need to change this for slab (-1/3  --> -1/2)
        end
        replace!(qPath, zeros(3) => zeros(3).+1e-7)
        if showProgress
                𝔻List = @showprogress 1. "Making Dynamical Matrices... " pmap(q -> 𝔻(q, crystal, couplings, charges, sumDepth, η), qPath)
        else
                𝔻List = pmap(q -> 𝔻(q, crystal, couplings, charges, sumDepth, η), qPath)
        end
        eigeninfo = eigen.(𝔻List)
        ω²Values = [round.(eig.values, digits=5) for eig in eigeninfo]
        fValues = map( x -> .√Complex.(x)./(2π), ω²Values)
        meVDispersion = 4.13567 .*fValues # convert THz to meV

        normalModes = [eig.vectors for eig in eigeninfo]

        return meVDispersion, normalModes
end


function getProjectedDispersion(qPath::Array, qzPath::Array, crystal::Crystal, couplings::Array)
        layerDispersions = []
        # Calculate dispersion for each layer to be projected onto surface BZ
        for qz in qzPath
                zLayerPath = map(q -> q + qz, qPath)
                zLayerDisp, = getDispersion(zLayerPath, crystal, couplings)
                push!(layerDispersions, zLayerDisp)
        end
        # Restructure output so that it can be plotted with plotDispersion(...)
        projectedDisp = Array{Array}(undef, length(qPath))
        for qi in eachindex(qPath)
            projectedDisp[qi] = []
            for layer in eachindex(layerDispersions)
                append!(projectedDisp[qi], layerDispersions[layer][qi])
            end
        end
        return projectedDisp
end


function getProjectedDispersion(qPath::Array, qzPath::Array, crystal::Crystal, couplings::Array, charges::Array, sumDepth::Int=5, η=nothing)
        if η == nothing
                η = 4*crystal.cellVol^(-1/3) # need to change this for slab (-1/3  --> -1/2)
        end
        replace!(qPath, zeros(3) => zeros(3).+1e-9)
        layerDispersions = []
        # Calculate dispersion for each layer to be projected onto surface BZ
        @showprogress 1 "Projecting dispersion... " for qz in qzPath
                zLayerPath = map(q -> q + qz, qPath)
                zLayerDisp, = getDispersion(zLayerPath, crystal, couplings, charges, sumDepth, η, false)
                push!(layerDispersions, zLayerDisp)
        end
        # Restructure output so that it can be plotted with plotDispersion(...)
        projectedDisp = Array{Array}(undef, length(qPath))
        for qi in eachindex(qPath)
            projectedDisp[qi] = []
            for layer in eachindex(layerDispersions)
                append!(projectedDisp[qi], layerDispersions[layer][qi])
            end
        end
        return projectedDisp
end



function plotDispersion(dispersion::Array, qPathParts::Array=[], qLabels::Array=[]; ylims::Array=[0.0, Inf], color::Symbol=:steelblue, lw::Real=1.5, size::Tuple=(600,350), title::String="")
    bands = []

    for i in 1:length(dispersion[1])
        band = [set[i] for set in dispersion]
        band = map(real, band)
        push!(bands, band)
    end
    plot(bands, xticks=(qPathParts, qLabels), legend=false, xtickfont=(13), size=size, color=color, lw=lw)
    if length(qPathParts) > 0
            plot!(qPathParts, seriestype=:vline, color=:black, linealpha=0.35)
            xlims!((1.0, qPathParts[end]))
    else
            xlims!((1.0, Inf))
    end
    ylims!(ylims[1], ylims[2])
    title!(title)
    ylabel!("Energy (meV)")
end
