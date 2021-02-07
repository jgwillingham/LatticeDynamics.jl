


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
