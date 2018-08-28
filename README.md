
[![Build Status](https://travis-ci.com/MurrellGroup/DPMeansClustering.jl.svg?branch=master)](https://travis-ci.com/MurrellGroup/DPMeansClustering.jl)

<a id='DPMeansClustering-1'></a>

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

Contains a variety of functions using dp based computation

## Installation
```julia
Pkg.clone("https://github.com/MurrellGroup/DPMeansClustering.jl.git")

```

## Set paths
```julia
using DPMeansClustering
```

## Run Tests
```julia
Pkg.test("DPMeansClustering")
