


using LinearAlgebra
using Printf

# Generalized dot product for any two equal size arrays
function dott(x::Array, y::Array)
        res = sum(x[i].*y[i] for i in eachindex(x))
        return res
end


function millerStringToArray(millerIndices::String)
        hkl = []
        sign = 1
        for char in millerIndices
                if char=="-"
                        sign = -1
                else
                        inx = parse(Int, char)
                        append!(hkl, sign*inx)
                        sign = 1
                end
        end
        return hkl
end


function getReciprocalVectors(latticeVectors::Array)
        a₁, a₂, a₃ = latticeVectors
        vol = abs(dot(a₁, cross(a₂, a₃))) # volume of unit cell
        b₁ = 2π/vol * cross(a₂, a₃)
        b₂ = 2π/vol * cross(a₃, a₁)
        b₃ = 2π/vol * cross(a₁, a₂)
        return [b₁, b₂, b₃]
end


function getSurfaceNormal(hkl::Array, reciprocalVectors::Array)
        G_hkl = dott(hkl, reciprocalVectors)
        surfaceNormal = G_hkl/norm(G_hkl)
        return surfaceNormal
end


function getAdaptedLatticeVectors(latticeVectors::Array, surfaceNormal::Array)
        searchRange = -2:2
        threshold = 10^-9
        a₁, a₂, a₃ = latticeVectors
        n = surfaceNormal
        # 1) within search range, collect lattice vectors that are parallel and non-parallel to the surface
        coplanars = []
        nonCoplanars = []
        for n₁ in searchRange
                for n₂ in searchRange
                        for n₃ in searchRange
                                Rℓ = n₁*a₁ + n₂*a₂ + n₃*a₃
                                if abs(dot(Rℓ, n)) < threshold && !(n₁==n₂==n₃==0)
                                        push!(coplanars, Rℓ)
                                end
                                if abs(dot(Rℓ, n) > threshold) && !(n₁==n₂==n₃==0)
                                        push!(nonCoplanars, Rℓ)
                                end
                        end
                end
        end
        # 2) sort to find the shortest in-plane lattice vectors
        coplanarLengthOrder = sortperm(map(norm, coplanars))
        orderedCoplanars = coplanars[coplanarLengthOrder]
        meshPrimitives = [ orderedCoplanars[1] ] # shortest in-plane lattice vector in new list 'meshPrimitives'
        # 3) sort to find the shortest out-of-plane lattice vector
        nonCoplanarLengthOrder = sortperm(map(norm, nonCoplanars))
        outOfPlaneVector = nonCoplanars[nonCoplanarLengthOrder[1]] # shortest out-of-plane lattice vector
        # 4) Find next shortest in-plane lattice vector which is independent of the other one
        vecAngle(v₁, v₂) = acos(dot(v₁, v₂) / (norm(v₁)*norm(v₂)) )
        isnonParallel(v₁, v₂) = ( mod(vecAngle(v₁, v₂), π) > 10^-5 )
        for Rℓ in orderedCoplanars
                if isnonParallel(Rℓ, meshPrimitives[1])
                        push!(meshPrimitives, Rℓ)
                        break
                end
        end
        adaptedLatticeVectors = [meshPrimitives[1], meshPrimitives[2], outOfPlaneVector]
        return adaptedLatticeVectors
end


function getSlabFractionalCoords(cartesianCoords::Array, adaptedLatticeVectors::Array, numCells::Int)
        a₁, a₂, a₃ = adaptedLatticeVectors
        vol = abs(dot(a₁, cross(a₂, a₃*numCells)))

        σ₁ = (1/vol) * cross(a₂, a₃*numCells)
        σ₂  = (1/vol) * cross(a₃*numCells, a₁)
        σ₃ = (1/vol) * cross(a₁, a₂)
        σ = [σ₁, σ₂, σ₃]
        fractionalCoords = [round(dot(cartesianCoords, σᵢ), digits=6) for σᵢ in σ]
        return fractionalCoords
end


function getSlabCell(bulkUnitCell::Array, latticeVectors::Array, adaptedLatticeVectors::Array, numCells::Int)
        outOfPlanePrimitive = adaptedLatticeVectors[3]
        numAtomsMovedToBottom = 0
        slabCell = []
        for ℓ3 in 1:numCells
                inx = @sprintf("_%0.0f", ℓ3) # index the atoms
                for atom in bulkUnitCell
                        element = atom[1]
                        bulkPosition_cartesian = dott(atom[2], latticeVectors)
                        position = bulkPosition_cartesian + (ℓ3-1)*outOfPlanePrimitive
                        fractionalCoords = getSlabFractionalCoords(position, adaptedLatticeVectors, numCells)
                        if (1.0 - abs(fractionalCoords[3])) < 10^-9
                                fractionalCoords[3] = 0.0
                                pushfirst!(slabCell, [element*"_0", fractionalCoords])
                                numAtomsMovedToBottom += 1
                        else
                                push!(slabCell, [element*inx, fractionalCoords])
                        end
                end
        end
        return slabCell, numAtomsMovedToBottom
end



function projectVector(vector::Vector{Real}, surfaceNormal::Vector{Real})
        n = surfaceNormal
        outOfPlane =  dot(vector, n) *  n
        inPlane = vector - outOfPlane
        return inPlane, outOfPlane
end





# Type for bulk crystal structural information
struct Crystal{T<:AbstractArray}
        unitCell::T
        latticeVectors::T
        reciprocalVectors::T

        function Crystal(unitCell, latticeVectors)
                reciprocalVectors = getReciprocalVectors(latticeVectors)
                new{AbstractArray}(unitCell, latticeVectors, reciprocalVectors)
        end
end


# Type for slab structural information
struct Slab{T<:AbstractArray}
        slabCell::T
        latticeVectors::T
        surface::String
        numCells::Int
        meshPrimitives::T
        meshReciprocals::T
        surfaceNormal::T
        numAtomsMovedToBottom::Int

        function Slab(bulkUnitCell, latticeVectors, surface, numCells)
                hkl = millerStringToArray(surface)
                reciprocalVectors = getReciprocalVectors(latticeVectors)
                surfaceNormal = getSurfaceNormal(hkl, reciprocalVectors)
                adaptedLatticeVectors = getAdaptedLatticeVectors(latticeVectors, surfaceNormal)
                meshPrimitives = adaptedLatticeVectors[1:2]
                outOfPlanePrimitive = adaptedLatticeVectors[3]
                meshReciprocals = getReciprocalVectors([meshPrimitives[1], meshPrimitives[2], surfaceNormal])
                slabCell, numAtomsMovedToBottom = getSlabCell(bulkUnitCell, latticeVectors, adaptedLatticeVectors, numCells)
                new{AbstractArray}(
                                slabCell,
                                adaptedLatticeVectors,
                                surface,
                                numCells,
                                meshPrimitives,
                                meshReciprocals,
                                surfaceNormal,
                                numAtomsMovedToBottom)
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
