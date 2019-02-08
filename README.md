# Detecting geometry via Betti curves

Repository for producing random sampling from unit balls in hyperbolic, Euclidian and spherical spaces, and analyzing them using tools from (CITE HERE).

The code is written in the Julia language using Julia v.0.6.4. 

`*.jl` files contain necessary functions (see descriptions below), and the notebooks illustrate the use and show the plots comparing the three geometries.

## Installation

This isn't a package (yet), so to use the functions, download the file into the working directory and run:
```
julia> include("[filename].jl")
```

## Sampling (`GeometricSampling.jl`)
`hyperbolic.jl` contains functions for sampling from different geometries and coomputing distances. Only sampling within a ball of radius 1 and in curvature 0, +1, -1 are implemented as of 02/08/2019.

### Sample ball
```
sample_ball(d::Int, numofpts=1, curvature=0.0, radius=1.0)::Array{Array{Float64,1},1}
```
Uniformly sample n=`numofpts` points in a ball of dimension `d` or radius `radius`, in a space with constant curvature `curvature`.

**Note:** This returns an array of coordinates of _the Euclidian representation_ of the points (i.e. points within the Poincare disk for `curvature` < 0, or points in the upper `d+1` hemisphere for `curvature` > 0).

```
sample_hyp(d::Int, numofpts=1, curvature=-1.0, radius=1.0)::Array{Array{Float64,1},1}
sample_euc(d::Int, numofpts=1, radius=1.0)::Array{Array{Float64,1},1}
sample_sph(d::Int, numofpts=1, curvature= 1.0, radius=1.0)::Array{Array{Float64,1},1}
```
Aliases for `sample_ball` in hyperbolic/Euclidian/spherical space respectively.

```
sample_sphere(d::Int, numofpts::Int=1, radius=1.0)::Array{Array{Float64,1},1}
```
Uniformly sample n=`numofpts` points on the surface of a sphere *in `d`-dimensional Euclidian space* (so a `(d-1)`-sphere)  of radius `radius`.

### Distances
```
distance_matrix(pts::Array{Array{Float64,1},1}; curvature=0.0)::Array{Float64, 2}
```
Returns the distance matrix for the collection of points `pts`, i.e. the `(i,j)`th entry is the distance between `pts[i]` and `pts[j]`, in the space with curvature `curvature`.

```
hyp_distance(pts::Array{Array{Float64,1},1}; curvature=-1.0)
euc_distance(pts::Array{Array{Float64,1},1})
sph_distance(pts::Array{Array{Float64,1},1}; curvature=1.0)
```
Aliases for `distance_matrix` in the hyperbolic, Euclidian, and spherical space respectively.

### Additional general-purpose functions

```
rejection_sampling(dens::Function, maxval::Float64, numofpts=1)::Array{Array{Float64,1},1}
```
Perform the simplest version of the rejection sampling of n=`numofpts` points from the density function `dens`, assuming that the support of the density function is `[0, maxval]`. The proposal distribution is taken to be the unifrom distribution on `[0, dens(maxval)]`.

```
to_density(matr::Array{T,2})
to_density!(matr::Array{T,2})
```
Convert a symmetric `n x n` matrix to a density matrix, i.e. replace `matr[i,j]` with the the number of entries in the upper triangle of `matr` that are less than `matr[i,j]`, divided by \binom{n}{2}. This ensures that the matrix entries are on `[0,1]` scale while preserving the Vietoris-Rips complex of the matrix.

## Betti curves
`AverageBettis.jl` contains functions for computing and plotting average Betti curves. We use Eirene to compute the Betti numbers.



