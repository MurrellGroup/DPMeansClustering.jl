
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

<a id='Functions-1'></a>
# Functions
**`DPMeansClustering.cluster`** &mdash; *Function*
```cluster(inputs, radius::Float64; distfunc=euclidean, center=mean, verbose::Int64=1,
            cycle_lim::Int64=30, triangle=false)
            ```
A wrapper that decides whether or not to use triangle inequality.

**`DPMeansClustering.get_initial_centroids`** &mdash; *Function*
```
get_initial_centroids(inputs, radius::Float64, input_length::Int64, distfunc, center, verbose::Int64)
```
Performs the first cycle of clustering in preparation for the creation of metacentroids

**`DPMeansClustering.get_initial_metacentroids`** &mdash; *Function*
```
get_initial_metacentroids(inputs, C_arr, MC_radius::Float64, distfunc, center, verbose::Int64)
```
Create metacentroids for the triangle inequality

**`DPMeansClustering.get_possible_centroids`** &mdash; *Function*
```
get_possible_centroids(MC_matrix, C_arr, R2MC_dist_arr::Array{Float64, 1},
                        C2MC_dist_matrix::Array{Array{Float64, 1}, 1}, radius::Float64, num_MC::Int64,
                        distfunc, MC_radius::Float64, verbose::Int64)
                        ```
Uses the triangle inequality in order to determine which centroids a read may cluster to. This is
done by finding the difference between the read to metacentroid(R2MC_dist_arr) as well as centroid
to metacentroid distance(C2MC_dist_matrix), since the triangle inequality mandates that the
distance from read to centroid must be at least the difference between the other two sides

Returns (possible_centroids_indices) which contains the indices of the centroids a read may belong
to. The indices correspond to the centroids' locations in the C_arr.         

**`DPMeansClustering.C2MC_dp_centers`** &mdash; *Function*
```
C2MC_dp_centers(inputs, radius::Float64; distfunc=euclidean, center=mean, verbose::Int64=1,
                    cycle_lim::Int64=30)
                    ```
                    Clusters the input kmer arrays using the triangle inequality for given radius


**`DPMeansClustering.dp_means`** &mdash; *Function*

```
dp_means(input_vectors::Array{Vector{Int}, 1}, radius::Float64; verbose = false)
```

Cluster `input_vectors` using euclidean distance metric and arithmetic mean, where a  vector with distance greater than `radius` from the nearest cluster forms a new cluster.

Returns (centroid_vectors, cluster_sizes, cluster_indices, centroid_indices) where cluster_indices is an array of arrays of indices of `input_vectors` grouped by cluster, and centroid_indices is an array of values where each value indexes the vector in `input_vectors` closest to the respective centroid.

Example:

```
inputs = [[5, 0, 0, 1, 1, 5],
          [5, 1, 0, 0, 1, 5],
          [5, 0, 1, 0, 1, 5],
          [0, 4, 6, 2, 0, 0],
          [0, 4, 6, 1, 1, 0],
          [0, 4, 6, 1, 0, 1]]
radius = 3.0
μs, sizes, indices, centroids = dp_means(inputs, radius)
μs == [[5, 1/3, 1/3, 1/3, 1, 5],
             [0, 4, 6, 4/3, 1/3, 1/3]]
sizes == [3, 3]
indices == [[1, 2, 3], [4, 5, 6]]
centroids == [1, 4]
```

**`DPMeansClustering.dp_centers`** &mdash; *Function*.



```
dp_centers(inputs, radius::Float64; distfunc = euclidean, center = mean, verbose = false, cycle_lim = 30)
```

Cluster `input_vectors` using given distance and mean calculations, where a  vector with distance greater than `radius` from the nearest cluster forms a new cluster. Runs a maximum of `cycle_lim` iterations.

Returns (centroid_vectors, cluster_sizes, cluster_indices, centroid_indices) where cluster_indices is an array of arrays of indices of `input_vectors` grouped by cluster, and centroid_indices is an array of values where each value indexes the vector in `input_vectors` closest to the respective centroid.

This is a generalization of `dp_means`.


```
dp_centers(inputs, radii::Vector{Float64}; distfunc = euclidean, center = mean, verbose = false, cycle_lim = [30...])
```

Cluster `input_vectors` using given distance and mean calculations recursively; `radii` should be an array of decreasing values, where each value is the radius of a successively finer clustering operation. During each clustering, a vector with distance greater than `radius` from the nearest cluster forms a new cluster. Runs a maximum of `cycle_lims` iterations during each cluster.

Returns (centroid_vectors, cluster_sizes, cluster_indices, centroid_indices) where cluster_indices is an array of arrays of indices of `input_vectors` grouped by cluster, and centroid_indices is an array of values where each value indexes the vector in `input_vectors` closest to the respective centroid.

This is a recursive implementation of `dp_centers(inputs, radius)` for faster clustering.
