


using LinearAlgebra: norm


function dott(x::Array, y::Array)
        res = sum(x[i].*y[i] for i in eachindex(x))
        return res
end


# Type for bulk crystal structural information
mutable struct Crystal{T<:AbstractArray}
        unitCell::T
        latticeVectors::T
        neighbors::Dict{String, AbstractArray}
        function Crystal(unitCell, latticeVectors)
                cartesian_unitCell = [[atom[1], dott(atom[2],latticeVectors)]
                                        for atom in unitCell]
                new{AbstractArray}(cartesian_unitCell, latticeVectors, Dict{String, AbstractArray}())
        end
end



function getNeighbors!(crystal::Crystal, threshold::Real, searchWidth::Integer=2)
        searchRange = -searchWidth:searchWidth
        neighbors = Dict{String, AbstractArray}()
        for atomᵢ in crystal.unitCell
                neighbors[atomᵢ[1]] = []
                rᵢ = atomᵢ[2]
                for atomⱼ in crystal.unitCell
                        xⱼ = atomⱼ[2]
                        for n₁ in searchRange
                                for n₂ in searchRange
                                        for n₃ in searchRange
                                                Rℓ = dott([n₁, n₂, n₃], crystal.latticeVectors)
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
        crystal.neighbors = neighbors
end
