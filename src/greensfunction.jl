


using LinearAlgebra


function blockSplit(matrix::Hermitian, blockSize::Int)
    matrixSize = size(matrix)[1]
    if mod(matrixSize, blockSize) != 0
        throw(ArgumentError("Block dimension must divide matrix dimension"))
    end
    blockIndices = [[(i:(i+blockSize-1), j:(j+blockSize-1)) for j in 1:blockSize:matrixSize] for i in 1:blockSize:matrixSize]
    # reduce(hcat, blockIndices) forms a matrix of indices such that matrix[inx...] has the same position in the matrix
    blockViews = map(x -> view(matrix, x...), reduce(hcat, blockIndices))
    return blockViews
end


function isTridiagonalSplit(blockViews::Array, tol::Real=10.0^-9.0)
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


function getPrincipalLayerSize(dynamicalMatrix::Hermitian, tol::Real=10.0^-9.0)
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
    throw(ArgumentError("Your dynamical matrix is not tridiagonal. Possibly it is too small."))
end


@inline function sanchoIterate(zI::Array, α::Array, β::Array, εˢ::Array, ε::Array)
    g = inv( zI - ε )
    newα = α*g*α
    newβ = β*g*β
    newεˢ = εˢ + α*g*β
    newϵ = ε + α*g*β + β*g*α
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
    # iterate
    counter = 0
    while counter <= iterNum
        α, β, εˢ, ε = sanchoIterate(zI, α, β, εˢ, ε)
        counter += 1
    # calculate surface Green's function
    Gω = inv( zI - εˢ )
    # Calculate spectral function / local density of states (LDOS) at (q,ω)
    Aω = (-1.0/π) * imag(tr(Gω))
    return Aω
end


function spectralFunction(qList::Array, εList::Array, crystal::Slab, couplings::Array; η::Real=10.0^-4, iterNum::Integer=15)
    εToω = 2π/4.13567
    ωList = εToω .* εList
end
