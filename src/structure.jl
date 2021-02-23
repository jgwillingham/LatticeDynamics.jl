


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
        threshold = 1e-9
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
        isnonParallel(v₁, v₂) = ( mod(vecAngle(v₁, v₂), π) > 1e-5 )
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
                        if (1.0 - abs(fractionalCoords[3])) < 1e-9
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


function getCartesianUnitCell(unitCell::Array, latticeVectors::Array)
        cartesianUnitCell = deepcopy(unitCell)
        for i in eachindex(unitCell)
                cartesianUnitCell[i][2] = dott(unitCell[i][2], latticeVectors)
        end
        return cartesianUnitCell
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
        cartesianUnitCell::T
        cellVol::Float64
        latticeVectors::T
        masses::T
        𝕄::T
        reciprocalVectors::T
        neighbors::Dict

        function Crystal(unitCell, latticeVectors, threshold)
                cartesianUnitCell = getCartesianUnitCell(unitCell, latticeVectors)
                cellVol = abs(dot(latticeVectors[1], cross(latticeVectors[2], latticeVectors[3])))
                masses = getMasses(unitCell)
                𝕄 = getMassMatrix(masses)
                reciprocalVectors = getReciprocalVectors(latticeVectors)
                neighbors = getBulkNeighbors(cartesianUnitCell, latticeVectors, threshold)
                new{AbstractArray}(unitCell,
                                cartesianUnitCell,
                                cellVol,
                                latticeVectors,
                                masses,
                                𝕄,
                                reciprocalVectors,
                                neighbors)
        end
end


# Type for slab structural information
struct Slab{T<:AbstractArray}
        unitCell::T
        cartesianUnitCell::T
        cellVol::Float64
        latticeVectors::T
        surface::String
        numCells::Int
        masses::T
        𝕄::T
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
                a₁, a₂, a₃ = adaptedLatticeVectors
                meshPrimitives = [a₁, a₂]
                cellVol = norm( cross(a₁, a₂) )
                outOfPlanePrimitive = a₃
                meshReciprocals = getReciprocalVectors([a₁, a₂, surfaceNormal])
                unitCell, numAtomsMovedToBottom = getSlabCell(bulkUnitCell, latticeVectors, adaptedLatticeVectors, numCells)
                cartesianUnitCell = getCartesianUnitCell(unitCell, [a₁, a₂, numCells*a₃] )
                neighbors = getSlabNeighbors(cartesianUnitCell, adaptedLatticeVectors, threshold, numCells)
                masses = getMasses(unitCell)
                𝕄 = getMassMatrix(masses)
                new{AbstractArray}(
                                unitCell,
                                cartesianUnitCell,
                                cellVol,
                                adaptedLatticeVectors,
                                surface,
                                numCells,
                                masses,
                                𝕄,
                                meshPrimitives,
                                meshReciprocals,
                                surfaceNormal,
                                neighbors,
                                numAtomsMovedToBottom)
        end
end




function getBulkNeighbors(cartesianUnitCell::Array, latticeVectors::Array, threshold::Real, searchWidth::Int=2)
        searchRange = -searchWidth:searchWidth
        a₁, a₂, a₃ = latticeVectors
        neighbors = Dict{String, AbstractArray}()
        for atomᵢ in cartesianUnitCell
                neighbors[atomᵢ[1]] = []
                rᵢ = atomᵢ[2]
                for atomⱼ in cartesianUnitCell
                        xⱼ = atomⱼ[2]
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


function getSlabNeighbors(cartesianUnitCell::Array, adaptedLatticeVectors::Array, threshold::Real, numCells::Int, searchWidth::Int=2)
        searchRange = -searchWidth:searchWidth
        a₁, a₂, a₃ = adaptedLatticeVectors
        slabPrimitives = [a₁, a₂, numCells*a₃]
        neighbors = Dict{String, AbstractArray}()
        for atomᵢ in cartesianUnitCell
                atomᵢLabel = atomᵢ[1]
                neighbors[atomᵢLabel] = []
                rᵢ = atomᵢ[2]
                for atomⱼ in cartesianUnitCell
                        atomⱼLabel = atomⱼ[1]
                        xⱼ = atomⱼ[2]
                        for n₁ in searchRange
                                for n₂ in searchRange
                                        Rℓ = n₁*a₁ + n₂*a₂
                                        rⱼ = Rℓ + xⱼ
                                        bondᵢⱼ = rⱼ - rᵢ
                                        bondLength = norm(bondᵢⱼ)
                                        if bondLength < threshold && bondLength > 0.0
                                                zFractionalCoord = getSlabFractionalCoords(rⱼ, adaptedLatticeVectors, numCells)[3]
                                                if abs(zFractionalCoord) < 1.0
                                                        push!(neighbors[atomᵢLabel], [atomⱼLabel, [bondᵢⱼ, Rℓ]])
                                                end
                                        end
                                end
                        end
                end
        end
        return neighbors
end



function getMassMatrix(masses::Array)
        Nₐ = 6.02e23 # Avagadro
        massList = [mass/Nₐ *1e21 for mass in masses] # atomic weights to single atom mass in 10^-24 kg = 10^-21 g
        diagBlocks = [repeat([1/√m], 3) for m in massList]
        diag = vcat(diagBlocks...)
        𝕄 = diagm(diag)
        return 𝕄
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
         masses = [atomicMassDict[split(element, "_")[1]] for element in atomsInUnitCell]
         return masses
end
