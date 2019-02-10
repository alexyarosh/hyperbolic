# Detecting geometry via Betti curves

Repository for producing random sampling from unit balls in hyperbolic, Euclidian and spherical spaces, and analyzing them using tools from (CITE HERE).

The code is written in the Julia language and requires Julia v.0.7. 

`*.jl` files contain necessary functions, and the notebooks illustrate the use and show the plots comparing the three geometries.

`GeometricSampling.jl` -- functions for sampling and distance 

`AverageBettis.jl` -- functions for computing and plotting average Betti curves

`sampling.ipynb` -- notebook illustrating sampling

`betti_curves.ipynb` -- notebook containing average Betti curves

## Installation

This isn't a package (yet), so to use the functions, download the file into the working directory and run:
```
julia> include("[filename].jl")
```

## Sampling (`GeometricSampling.jl`)
`GeometricSampling.jl` contains functions for sampling from different geometries and coomputing distances. Only sampling within a ball of radius 1 and in curvature 0, +1, -1 are implemented as of 02/08/2019.

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

## Betti curves (`AverageBettis.jl`)
`AverageBettis.jl` contains functions for computing and plotting average Betti curves. 

```
bettis(matr::Array{Float64,2}, 
       maxdim::Int; 
       mintime::Float64 = -Inf,
       maxtime::Float64 = Inf
       numofsteps::Int = Inf
       method::Symbol = :ripser)::Array{Float64,2}
```       
Compute Betti numbers of the Vietoris-Rips complex defined by the distance matrix `matr` up to dimension `maxdim`. Parameters `mintime`, `maxtime`, `numofsteps` define the filtration, where `mintime` and `maxtime` is the start and end points of the filtration, and `numofsteps` is the number of steps in the filtration. Parameter `method` determines the software to use for persistent homology (possible values `method=:ripser`  for `Ripser.jl` and `method=:eirene` for `Eirene.jl`).

Returns  `numofsteps x maxdim` matrix, where `(i,j)`'th entry is the `j`th Betti number at filtration step `i`. Betti_0 is discarded.

```
average_bettis(arrs::Array{Array{Float64,2},1})::Array{Float64,2}
```
Assuming `arrs` is array of outputs of `bettis(..)` of same size, find average values for Betti numbers for each point in the filtration for each dimesion.

```
std_bettis(arrs::Array{Array{Float64,2},1})::Array{Float64,2}
```
Assuming `arrs` is array of outputs of `bettis(..)` of same size, find standard deviations for Betti numbers for each point in the filtration for each dimesion.

```
plot_averages(xvals::Array{Float64,1}, 
              means::Array{Float64,1}, 
              stds::Array{Float64,1}; 
              ribbon::Bool=true, 
              label::String = "", 
              linestyle = :solid, 
              color = :auto)
              
plot_averages!(xvals, means, stds; ribbon=true, label="", linestyle=:solid, color=:auto)
```
Plot average Betti curves `means` with standard deviations `stds` at filtration values given by `xvals`. Ribbon parameter determins whether to plot the error ribbon (of width=standard deviation) around the curve.


```
plot_averages(xvals::Array{Float64,1}, 
              means::Array{Float64,1}, 
              stds::Array{Float64,1}; 
              ribbon::Bool=true, 
              label::String = "", 
              linestyle = :solid, 
              color = :auto)
              
plot_averages!(xvals, means, stds; ribbon=true, label="", linestyle=:solid, color=:auto)
```
Plot average Betti curves `means` with standard deviations `stds` at filtration values given by `xvals`. Ribbon parameter determins whether to plot the error ribbon (of width=standard deviation) around the curve.

```
plot_averages(xvals::Array{Float64,1}, 
              file::String; 
              dim=1, 
              ribbon=true, 
              label="", 
              linestyle=:solid, 
              color=:auto)

plot_averages!(xvals, file::String; dim=1, ribbon=true, label="", linestyle=:solid, color=:auto)
```
Plot average Betti curves in dimension `dim` at filtration values given by `xvals`, given that `file` contains _all_ the Betti numbers (e.g. it contains many outputs of `bettis`).
