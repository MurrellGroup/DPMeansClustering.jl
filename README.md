# DPMeansClustering

- [`DPMeansClustering.cluster`]
- [`DPMeansClustering.get_initial_centroids`]
- [`DPMeansClustering.get_initial_metacentroids`]
- [`DPMeansClustering.get_possible_centroids`]
- [`DPMeansClustering.C2MC_dp_centers`]
- [`DPMeansClustering.dp_centers`]
- [`DPMeansClustering.dp_centers`]
- [`DPMeansClustering.dp_means`]

## DOCS
https://murrellgroup.github.io/DPMeansClustering.jl/

## Synopsis

Like K-means, but instead of having to choose the number of clusters (which you usually don't know), you pick a cluster radius. Points falling outside the radius of any current clusters spawn new clusters.

## Installation
```julia
using Pkg
Pkg.add(url="https://github.com/MurrellGroup/DPMeansClustering.jl")
```
