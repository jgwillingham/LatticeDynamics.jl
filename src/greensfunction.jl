

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


function isTridiagonalSplit(blockViews::Array, tol::Float64=1e-9)
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


function getPrincipalLayerSize(dynamicalMatrix::Hermitian, tol::Float64=1e-9)
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


@inline function sanchoIterate!(α::T, β::T, εˢ::T, ε::T, zI::T) where T<:Union{Array, SubArray}
    g = inv(zI - ε)
    εˢ[1:end,1:end] = εˢ + α*g*β
    ε[1:end,1:end] = ε + α*g*β + β*g*α
    α[1:end,1:end] = α*g*α
    β[1:end,1:end] = β*g*β
end


function getLDOS(ω::Float64, η::Float64, Dblocks::Array, iterNum::Int)
    # initial principal layer blocks from dynamical matrix
    α = Matrix(Dblocks[1,2])
    β = Matrix(Dblocks[2,1])
    εˢ = Matrix(Dblocks[1,1])
    ε = Matrix(Dblocks[2,2])

    z = ω^2 + im*η
    zI = z*Matrix(I, size(ε))
    # iterate / decimate
    counter = 1
    while counter <= iterNum
        sanchoIterate!(α, β, εˢ, ε, zI)
        counter += 1
    end
    # surface Green's function:
    Gω = inv(zI - εˢ)
    # bulk Green's function: (yields bulk projection without surface states)
    #Gω = inv(zI - ε)

    Aω = (-1.0/π) * imag(tr(Gω))

    return Aω
end


function getSpectrum(qList::Array, εList::Array{Float64,1}, crystal::Slab, couplings::Array; imag::Float64=1e-4, iterNum::Int=22)
    εToω = 2π/4.13567
    ωList = εToω .* εList

    testq = [0.1, 0.1, 0.1]
    testD = 𝔻(testq, crystal, couplings)
    PLSize = getPrincipalLayerSize(testD)
    atomDepth = 2*PLSize ÷ 3 # minimum number of atoms to consider to get necessary D blocks

    Aqω = Array{Array,1}(undef, length(qList))
    DList = @showprogress 0 "Building Dynamical Matrices... " pmap(x -> 𝔻(x, crystal, couplings, atomDepth), qList)
    DBlockList = map(x-> blockSplit(x, PLSize), DList)
    @showprogress 1 "Calculating LDOS... " for i in eachindex(DBlockList)
        energyCurve = pmap(x -> getLDOS(x, imag, DBlockList[i], iterNum), ωList)
        Aqω[i] = energyCurve
    end
    return Aqω
end


function getSpectrum(qList::Array, εList::Array{Float64,1}, crystal::Slab, couplings::Array, charges::Array, sumDepth::Int=5, η=nothing; imag::Float64=1e-4, iterNum::Int=22)
    εToω = 2π/4.13567
    ωList = εToω .* εList

    if η == nothing
            η = 4*crystal.cellVol^(-1/3) # need to change this for slab (-1/3  --> -1/2)
    end
    replace!(qList, zeros(3) => zeros(3).+1e-9)

    testq = rand(3)
    testD = 𝔻(testq, crystal, couplings, charges, sumDepth, η)
    PLSize = getPrincipalLayerSize(testD)
    println("Principal layer: ", PLSize ÷ 3, " atoms" )
    atomDepth = 2*PLSize ÷ 3 # minimum number of atoms to consider to get necessary D blocks

    Aqω = Array{Array,1}(undef, length(qList))
    DList = @showprogress 0 "Building Dynamical Matrices... " pmap(x -> 𝔻(x, crystal, couplings, charges, sumDepth, η, atomDepth), qList)
    DBlockList = map(x-> blockSplit(x, PLSize), DList)
    @showprogress 1 "Calculating LDOS... " for i in eachindex(DBlockList)
        energyCurve = pmap(x -> getLDOS(x, imag, DBlockList[i], iterNum), ωList)
        Aqω[i] = energyCurve
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


function getEnergySurface(energy::Float64, qPath1::Array, qPath2::Array, crystal::Slab, couplings::Array; η::Float64=1e-4, iterNum::Int=22)
    εToω = 2π/4.13567
    ω = εToω * energy

    testq = [0.1, 0.2, 0.3]
    testD = 𝔻(testq, crystal, couplings)
    PLSize = getPrincipalLayerSize(testD)
    atomDepth = 2*PLSize ÷ 3

    qArray = [[q₁+q₂ for q₂ in qPath2] for q₁ in qPath1]
    DArray = @showprogress 1 "Building Dynamical Matrices... " map(x -> 𝔻(x, crystal, couplings, atomDepth=atomDepth), reduce(hcat, qArray) )
    DBlocksArray = map(x->blockSplit(x, PLSize), DArray)
    εLDOS = @showprogress 1 "Calculating LDOS... " map(x->getLDOS(ω,η,x,iterNum), DBlocksArray)
    return εLDOS
end


function plotEnergySurface(εLDOS::Array)
    heatmap(log10.(εLDOS), xticks=([],[]), yticks=([],[]))
end
