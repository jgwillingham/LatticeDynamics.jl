


using LinearAlgebra


function blockSplit(matrix::Matrix, blockSize::Int)
    matrixSize = size(Matrix)[1]
    if mod(matrixSize, blockSize) != 0
        throw(ArgumentError("Block dimension must divide matrix dimension"))
    end
    blockIndices = [(i:i+2, j:j+2) for i in 1:3:matrixSize for j in 1:3:matrixSize]
    blockViews = [@view matrix[inx...] for inx in blockIndices]
    return blockViews
end


function isTridiagonal(blockViews::Array)
    absAverages = []
end


function getPrincipalLayerSize(dynamicalMatrix::Matrix)

end
