


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



function projectVector(vector::Array, surfaceNormal::Array)
        n = surfaceNormal
        outOfPlane =  dot(vector, n) *  n
        inPlane = vector - outOfPlane
        return inPlane, outOfPlane
end





# Type for bulk crystal structural information
struct Crystal{T<:AbstractArray}
        unitCell::T
        latticeVectors::T
        masses::T
        reciprocalVectors::T
        neighbors::Dict

        function Crystal(unitCell, latticeVectors, threshold)
                masses = getMasses(unitCell)
                reciprocalVectors = getReciprocalVectors(latticeVectors)
                neighbors = getBulkNeighbors(unitCell, latticeVectors, threshold)
                new{AbstractArray}(unitCell, latticeVectors, masses, reciprocalVectors, neighbors)
        end
end


# Type for slab structural information
struct Slab{T<:AbstractArray}
        unitCell::T
        latticeVectors::T
        surface::String
        numCells::Int
        masses::T
        meshPrimitives::T
        meshReciprocals::T
        surfaceNormal::T
        neighbors::Dict
        numAtomsMovedToBottom::Int

        function Slab(bulkUnitCell, latticeVectors, surface, numCells, threshold)
                hkl = millerStringToArray(surface)
                reciprocalVectors = getReciprocalVectors(latticeVectors)
                surfaceNormal = getSurfaceNormal(hkl, reciprocalVectors)
                adaptedLatticeVectors = getAdaptedLatticeVectors(latticeVectors, surfaceNormal)
                meshPrimitives = adaptedLatticeVectors[1:2]
                outOfPlanePrimitive = adaptedLatticeVectors[3]
                meshReciprocals = getReciprocalVectors([meshPrimitives[1], meshPrimitives[2], surfaceNormal])
                unitCell, numAtomsMovedToBottom = getSlabCell(bulkUnitCell, latticeVectors, adaptedLatticeVectors, numCells)
                neighbors = getSlabNeighbors(unitCell, adaptedLatticeVectors, threshold)
                masses = getMasses(unitCell)
                new{AbstractArray}(
                                unitCell,
                                adaptedLatticeVectors,
                                surface,
                                numCells,
                                masses,
                                meshPrimitives,
                                meshReciprocals,
                                surfaceNormal,
                                neighbors,
                                numAtomsMovedToBottom)
        end
end




function getBulkNeighbors(unitCell::Array, latticeVectors::Array, threshold::Real, searchWidth::Integer=2)
        searchRange = -searchWidth:searchWidth
        a₁, a₂, a₃ = latticeVectors
        neighbors = Dict{String, AbstractArray}()
        for atomᵢ in unitCell
                neighbors[atomᵢ[1]] = []
                rᵢ = dott(atomᵢ[2], latticeVectors)
                for atomⱼ in unitCell
                        xⱼ = dott(atomⱼ[2], latticeVectors)
                        for n₁ in searchRange
                                for n₂ in searchRange
                                        for n₃ in searchRange
                                                Rℓ = n₁*a₁ + n₂*a₂ + n₃*a₃
                                                rⱼ  = Rℓ + xⱼ
                                                bondᵢⱼ = rⱼ - rᵢ
                                                bondLength = norm(bondᵢⱼ)
                                                if bondLength < threshold && bondLength !=0.0
                                                        push!(neighbors[atomᵢ[1]], [atomⱼ[1], [bondᵢⱼ, Rℓ]])
                                                end
                                        end
                                end
                        end
                end
        end
        return neighbors
end


function getSlabNeighbors(unitCell::Array, adaptedLatticeVectors::Array, threshold::Real, searchWidth::Integer=2)
        searchRange = -searchWidth:searchWidth
        a₁, a₂, a₃ = adaptedLatticeVectors
        neighbors = Dict{String, AbstractArray}()
        for atomᵢ in unitCell
                atomᵢLabel = atomᵢ[1]
                neighbors[atomᵢLabel] = []
                rᵢ = dott(atomᵢ[2], adaptedLatticeVectors)
                for atomⱼ in unitCell
                        atomⱼLabel = atomⱼ[1]
                        xⱼ = dott(atomⱼ[2], adaptedLatticeVectors)
                        for n₁ in searchRange
                                for n₂ in searchRange
                                        Rℓ = n₁*a₁ + n₂*a₂
                                        rⱼ = Rℓ + xⱼ
                                        bondᵢⱼ = rⱼ - rᵢ
                                        bondLength = norm(bondᵢⱼ)
                                        zFractionalCoord = atomⱼ[2][3]
                                        if 10^-9 < bondLength < threshold && 0 <= abs(zFractionalCoord) < 1
                                                push!(neighbors[atomᵢLabel], [atomⱼLabel, [bondᵢⱼ, Rℓ]])
                                        end
                                end
                        end
                end
        end
        return neighbors
end


function getMasses(unitCell::Array)
        atomsInUnitCell = [atom[1] for atom in unitCell]
        atomicMassDict = Dict("H"=> 1.008,
                         "He"=> 4.002602,
                         "Li"=> 6.94,
                         "Be"=> 9.0121831,
                         "B"=> 10.81,
                         "C"=> 12.011,
                         "N"=> 14.007,
                         "O"=> 15.999,
                         "F"=> 18.998403163,
                         "Ne"=> 20.1797,
                         "Na"=> 22.98976928,
                         "Mg"=> 24.305,
                         "Al"=> 26.9815385,
                         "Si"=> 28.085,
                         "P"=> 30.973761998,
                         "S"=> 32.06,
                         "Cl"=> 35.45,
                         "Ar"=> 39.948,
                         "K"=> 39.0983,
                         "Ca"=> 40.078,
                         "Sc"=> 44.955908,
                         "Ti"=> 47.867,
                         "V"=> 50.9415,
                         "Cr"=> 51.9961,
                         "Mn"=> 54.938044,
                         "Fe"=> 55.845,
                         "Co"=> 58.933194,
                         "Ni"=> 58.6934,
                         "Cu"=> 63.546,
                         "Zn"=> 65.38,
                         "Ga"=> 69.723,
                         "Ge"=> 72.63,
                         "As"=> 74.921595,
                         "Se"=> 78.971,
                         "Br"=> 79.904,
                         "Kr"=> 83.798,
                         "Rb"=> 85.4678,
                         "Sr"=> 87.62,
                         "Y"=> 88.90584,
                         "Zr"=> 91.224,
                         "Nb"=> 92.90637,
                         "Mo"=> 95.95,
                         "Tc"=> 97.90721,
                         "Ru"=> 101.07,
                         "Rh"=> 102.9055,
                         "Pd"=> 106.42,
                         "Ag"=> 107.8682,
                         "Cd"=> 112.414,
                         "In"=> 114.818,
                         "Sn"=> 118.71,
                         "Sb"=> 121.76,
                         "Te"=> 127.6,
                         "I"=> 126.90447,
                         "Xe"=> 131.293,
                         "Cs"=> 132.90545196,
                         "Ba"=> 137.327,
                         "La"=> 138.90547,
                         "Ce"=> 140.116,
                         "Pr"=> 140.90766,
                         "Nd"=> 144.242,
                         "Pm"=> 144.91276,
                         "Sm"=> 150.36,
                         "Eu"=> 151.964,
                         "Gd"=> 157.25,
                         "Tb"=> 158.92535,
                         "Dy"=> 162.5,
                         "Ho"=> 164.93033,
                         "Er"=> 167.259,
                         "Tm"=> 168.93422,
                         "Yb"=> 173.045,
                         "Lu"=> 174.9668,
                         "Hf"=> 178.49,
                         "Ta"=> 180.94788,
                         "W"=> 183.84,
                         "Re"=> 186.207,
                         "Os"=> 190.23,
                         "Ir"=> 192.217,
                         "Pt"=> 195.084,
                         "Au"=> 196.966569,
                         "Hg"=> 200.592,
                         "Tl"=> 204.38,
                         "Pb"=> 207.2,
                         "Bi"=> 208.9804,
                         "Po"=> 209.0,
                         "At"=> 210.0,
                         "Rn"=> 222.0,
                         "Fr"=> 223.0,
                         "Ra"=> 226.0,
                         "Ac"=> 227.0,
                         "Th"=> 232.0377,
                         "Pa"=> 231.03588,
                         "U"=> 238.02891,
                         "Np"=> 237.0,
                         "Pu"=> 244.0,
                         "Am"=> 243.0,
                         "Cm"=> 247.0,
                         "Bk"=> 247.0,
                         "Cf"=> 251.0,
                         "Es"=> 252.0,
                         "Fm"=> 257.0,
                         "Md"=> 258.0,
                         "No"=> 259.0,
                         "Lr"=> 262.0,
                         "Rf"=> 267.0,
                         "Db"=> 268.0,
                         "Sg"=> 271.0,
                         "Bh"=> 274.0,
                         "Hs"=> 269.0,
                         "Mt"=> 276.0,
                         "Ds"=> 281.0,
                         "Rg"=> 281.0,
                         "Cn"=> 285.0,
                         "Nh"=> 286.0,
                         "Fl"=> 289.0,
                         "Mc"=> 288.0,
                         "Lv"=> 293.0,
                         "Ts"=> 294.0,
                         "Og"=> 294.0)
         masses = [atomicMassDict[element] for element in atomsInUnitCell]
         return masses
end
