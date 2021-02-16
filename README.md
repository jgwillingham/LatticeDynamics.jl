# LatticeDynamics.jl

This package provides basic functions for modeling harmonic lattice dynamics in bulk crystals as well as in slab or semi-infinite geometries.

To construct a model, a user inputs two pieces of information:

* Basic crystal information (*e.g. lattice vectors, unit cell, surface miller indices, etc.*)
* A parameter set which determines the interatomic force constants

With this, the dynamical matrix can be constructed so that the phonon spectrum and normal modes are readily calculated. The motivation for making the package comes from studying topological phonons so there is a mild bias toward functionality that lends itself to that field. The goal is to be able to seemlessly explore bulk and surface modes of a model with clean succinct code. 
_______________________________
#### Contents
* [Crystal Structure](#Crystal-Structure)
    - [Bulk Crystals](#Bulk-Crystals)
    - [Surfaces](#Surfaces)
* [Interactions](#Interactions)
    - [Short-Range Forces](#Short-Range-Forces)
    - [Coulomb Effective Charges](#Coulomb-Effective-Charges)
* [The Phonon Spectrum](#The-Phonon-Spectrum)
    - [Paths Through the (Surface) BZ](#Paths-Through-the-Surface-BZ)
    - [Calculating/Plotting Dispersion](#CalculatingPlotting-Dispersion)
    - [Green's Function Methods](#Greens-Function-Methods)
________________________________
## Crystal Structure

### Bulk Crystals
For a bulk crystal, we need 
1) Lattice vectors
2) Information about the unit cell
3) A radius which captures all interacting neighbors

Consider the following schematic example:
```julia
using LatticeDynamics

# For lattice vectors a1, a2, a3
latticevectors = [a1, a2, a3]

# For two atoms in the unit cell:
atom2info = ["Element1", frac_coords1] # "Element1" is the element symbol, e.g. "Al", "C", etc.
atom2info = ["Element2", frac_coords2]
unitcell = [atom1info, atom2info]

# Neighbor threshold in angstroms  <- maximum extent of the short-range interactions
threshold = 5.0

mycrystal = Crystal(unitcell, latticevectors, threshold)
```
Now `mycrystal` contains all the structural information about the crystal including reciprocal lattice information, atomic masses, and nearest neighbors for each atom.

### Surfaces

To study phonons at the surface of a crystal, we need a `Slab`. To make a slab, we use the same 3 pieces of bulk crystal information along with

4) Surface miller indices
5) The thickness of the slab. 

Here is an example:
```julia
surface = "hkl" # This could be "110", "001", etc.
thickness = 40 # Number of bulk unit cells to stack 
myslab = Slab(unitcell, latticevectors, surface, thickness, threshold)
```
**Note that the miller indices `"hkl"` must be with respect to the reciprocal lattice vectors derived from the given lattice vectors.** This can lead to potentially unexpected results. Recall that a surface labeled by "hkl" means that the vector `G=h*b1 + k*b2 + l*b3` is normal to the surface where `b1`, `b2` and `b3` are reciprocal lattice vectors derived from the direct lattice vectors in the usual way. Care must be taken to ensure the miller indices `"hkl"` given to `Slab` match the expected reciprocal lattice vectors. One can always look at `myslab.surfaceNormal` to ensure it is pointed in the proper direction.

_______________________________

## Interactions
Right now, LatticeDynamics.jl supports two types of interactions:
* Short-range radial interactions
* Long-range Coulomb interactions

### Short-Range Forces

The short-range forces are modeled by a radial potential. So the interaction between atoms at positions ***r***<sub>i</sub> and ***r***<sub>j</sub> looks like *V*(|***r***<sub>i</sub> - ***r***<sub>j</sub>|). With this assumption, the force constant matrix only depends on the value of the first and second derivatives at the equilibrium separation. These two values are taken to be phenomenological parameters which we might call *A*<sub>ij</sub> and *B*<sub>ij</sub>. Packaging these as a tuple (*A*<sub>ij</sub>, *B*<sub>ij</sub>), we organize the short-range interaction parameters for every pair of atom types as 

```julia
# For two atoms in the unit cell:
shortrange_couplings = [[ (A11, B11) , (A12, B12) ],
                        [ (A12, B12) , (A22, B22) ]]
```
So the short-range interaction parameters between atoms *i* and *j* is just obtained by `shortrange_couplings[i][j]`. To translate this to the couplings for a similar slab, we use the `getSlabCouplingArray` function
```julia
slabcouplings = getSlabCouplingArray(myslab, shortrange_couplings)
```

### Coulomb Effective Charges

For the long-range Coulomb interactions, the form of the potential is of course 1/r. So we simply assign to each atom in the unit cell an effective charge *Z*<sub>i</sub> which is understood to be in units of the fundamental electron charge *e*. So we write
```julia
# Again for two atoms in the unit cell:
charges = [Z1, Z2]
@assert sum(charges) == 0.0
```
Note that the net charge in the unit cell must be zero.

____________________________________

## The Phonon Spectrum

### Paths Through the (Surface) BZ

We can now calculate the phonon spectrum. Given the high symmetry points in the Brillouin zone, we can build a path connecting them with the `buildPath` function. For example, if we have high symmetry points Γ, X, W, and K, we can build the path Γ -> X -> W -> K -> Γ, but we have to specify how many points to include. We do this with a `pointdensity`. We then build the path with 
```julia
pointdensity = 35 # sampling rate <- how many points per angstrom should we calculate the energies
bzpath, bzpathparts = buildPath([Γ, X, W, K, Γ], pointdensity)
```
When studying surfaces it is often convenient to use the `projectVector` function which can project a high symmetry point of the 3D BZ to one in the surface BZ
```julia
inplaneX, outofplaneX = projectVector(X, myslab.surfaceNormal)
```

### Calculating/Plotting Dispersion

With the Brillouin zone path defined, we can get the dispersion with `getDispersion` and plot it with `plotDispersion`
```julia
disp = getDispersion(bzpath, mycrystal, couplings)
labels = ["Γ", "X", "W", "K", "Γ"]
plotDispersion(disp, bzpathparts, labels)
```
This calculates the phonon dispersion along the given path Γ -> X -> W -> K -> Γ and plots it. The `bzpathparts` is just for plotting. It tells where on the plot to show the high symmetry points. 

### Green's Function Methods

 The `getDispersion` function can be used to directly diagonalize a slab dynamical matrix, but it is often useful to consider an alternative method to obtaining the surface spectrum; namely the recursion-decimation Green's function algorithm (Sancho 1985). This method has the advantage that in its output, it is easier to distinguish between bulk-projected modes and surface resonances. It is called with the `getSpectrum` function. Here is an example:
 ```julia
 # construct BZ path and list of energies to calculate the LDOS for
 sbzpath, sbzpathparts = buildPath(highsymmpoints, pointdensity)
 energylist = range(0., stop=20., length=750)
 
 # calculate the spectral function for the given wavevectors and energies
 spec = getSpectrum(sbzpath, energylist, myslab, slabcouplings)
 plotSpectrum(spec, energylist, sbzpathparts)
 ```
 Note that even though a slab is given to this function, it uses information so that it effectively considers a semi-infinite crystal. 
 
 Also if there is a single energy of interest (for example the energy of a Weyl point), a plot of the isoenergy surface over some region of the surface BZ can be made with the `getEnergySurface` and `plotEnergySurface` functions.
