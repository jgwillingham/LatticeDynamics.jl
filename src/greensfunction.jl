


using LinearAlgebra
using ProgressMeter


function blockSplit(matrix::Hermitian, blockSize::Int)
    matrixSize = size(matrix)[1]
    if mod(matrixSize, blockSize) != 0
        throw(ArgumentError("Block dimension must divide matrix dimension"))
    end
    blockIndices = [[(i:(i+blockSize-1), j:(j+blockSize-1)) for i in 1:blockSize:matrixSize] for j in 1:blockSize:matrixSize]
    # reduce(hcat, blockIndices) forms a matrix of indices such that matrix[inx...] has the same position in the matrix
    blockViews = map(x -> view(matrix, x...), reduce(hcat, blockIndices))
    return blockViews
end


function isTridiagonalSplit(blockViews::Array, tol::Real=1e-9)
    absBlocks = map(x->abs.(x), blockViews)
    # calculate the average abs value of elements in each 3x3 block
    absAverages = map(x -> sum(x)/9.0, absBlocks)
    # check if all blocks off of the tridiagonal spine are sufficiently close to zero
    # this takes advantage of the matrix being Hermitian by checking only the upper triangle
    for i in 1:size(absAverages)[1]
        if all(x -> x < tol, absAverages[i, i+2:end])
            continue
        else
            return false
        end
    end
    return true
end


function getPrincipalLayerSize(dynamicalMatrix::Hermitian, tol::Real=1e-9)
    Dsize = size(dynamicalMatrix)[1]
    maxBlockSize = Dsize ÷ 3 # principal layers this big will involve interacting surfaces
    for blockSize in 3:3:maxBlockSize
        if mod(Dsize, blockSize) == 0
            DblockViews = blockSplit(dynamicalMatrix, blockSize)
            if isTridiagonalSplit(DblockViews, tol)
                return blockSize
            end
        end
    end
    throw(ArgumentError("Your dynamical matrix is either too small or not block-tridiagonal"))
end


@inline function sanchoIterate(zI::Array, α::AbstractArray, β::AbstractArray, εˢ::AbstractArray, ε::AbstractArray)
    g = inv(zI - ε)
    newα = α*g*α
    newβ = β*g*β
    newεˢ = εˢ + α*g*β
    newε = ε + α*g*β + β*g*α
    return newα, newβ, newεˢ, newε
end


function getLDOS(ω::Real, η::Real, Dblocks::Array, iterNum::Integer)
    # initial principal layer blocks from dynamical matrix
    α = Dblocks[1,2]
    β = Dblocks[2,1]
    εˢ = Dblocks[1,1]
    ε = Dblocks[2,2]

    z = ω^2 + im*η
    blockSize = size(εˢ)[1]
    zI = z*Matrix(I, blockSize, blockSize)
    # iterate / decimate
    counter = 1
    while counter <= iterNum
        α, β, εˢ, ε = sanchoIterate(zI, α, β, εˢ, ε)
        counter += 1
    end
    # surface Green's function:
    Gω = inv(zI - εˢ)
    # bulk Green's function: (yields bulk projection without surface states)
    #Gω = inv(zI - ε)

    Aω = (-1.0/π) * imag(tr(Gω))

    return Aω
end


function getSpectrum(qList::Array, εList::Array, crystal::Slab, couplings::Array; η::Real=1e-4, iterNum::Integer=22)
    εToω = 2π/4.13567
    ωList = εToω .* εList

    testq = [0.1, 0.1, 0.1]
    testD = 𝔻(testq, crystal, couplings)
    PLSize = getPrincipalLayerSize(testD)

    Aqω = []
    prog = Progress(length(qList), 1, "", 50) # makes a progress bar
    for q in qList
        dynamicalMatrix = 𝔻(q, crystal, couplings)
        blocks = blockSplit(dynamicalMatrix, PLSize)
        energyCurve = map(x -> getLDOS(x, η, blocks, iterNum), ωList)
        push!(Aqω, energyCurve)
        next!(prog)
    end
    return Aqω
end


function plotSpectrum(Aqω::Array, εList::Array=[], qPathParts::Array=[], qLabels::Array=[]; size::Tuple=(600, 350), numyticks::Integer=5, color::Symbol=:inferno, title::String="")
    N = numyticks
    if length(εList)!=0
        yticklocs = [n/N*length(εList) for n in 0:N]
        ytickvals = [round(n/N*max(εList...), digits=2) for n in 0:N]
        yticks = (yticklocs, ytickvals)
    else
        yticks=([],[])
    end

    plottableAqω = log10.(hcat(Aqω...))
    heatmap(plottableAqω, xticks=(qPathParts, qLabels), xtickfont=(13), yticks=yticks, size=size, color=color, title=title, ylabel="Energy (meV)")
    if length(qPathParts) > 0
        plot!(qPathParts, seriestype=:vline, color=:white, linealpha=0.35, legend=false)
    end
end


function getEnergySurface(energy::Float64, qPath1::Array, qPath2::Array, crystal::Slab, couplings::Array; η::Float64=1e-4, iterNum::Integer=22)
    εToω = 2π/4.13567
    ω = εToω * energy

    testq = [0.1, 0.1, 0.1]
    testD = 𝔻(testq, crystal, couplings)
    PLSize = getPrincipalLayerSize(testD)

    qArray = [[q₁+q₂ for q₂ in qPath2] for q₁ in qPath1]
    DArray = map(x -> 𝔻(x, crystal, couplings), reduce(hcat, qArray) )
    DBlocksArray = map(x->blockSplit(x, PLSize), DArray)
    εLDOS = map(x->getLDOS(ω,η,x,iterNum), DBlocksArray)
    return εLDOS
end


function plotEnergySurface(εLDOS::Array)
    heatmap(log10.(εLDOS), xticks=([],[]), yticks=([],[]))
end
