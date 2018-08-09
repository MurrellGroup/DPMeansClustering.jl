
<a id='Clustering-Functions-1'></a>

# Clustering Functions

- [`DPMeansClustering.cluster`](dp.md#DPMeansClustering.cluster)
- [`DPMeansClustering.dp_centers`](dp.md#DPMeansClustering.dp_centers)
- [`DPMeansClustering.dp_means`](dp.md#DPMeansClustering.dp_means)

<a id='DPMeansClustering.dp_means' href='#DPMeansClustering.dp_means'>#</a>
**`DPMeansClustering.dp_means`** &mdash; *Function*.



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


<a target='_blank' href='https://github.com/MurrellGroup/DPMeansClustering.jl/blob/25d48980e173119cfda767725fe058da77c60d99/src/dp.jl#L297-L325' class='documenter-source'>source</a><br>

<a id='DPMeansClustering.dp_centers' href='#DPMeansClustering.dp_centers'>#</a>
**`DPMeansClustering.dp_centers`** &mdash; *Function*.



```
dp_centers(inputs, radius::Float64; distfunc = euclidean, center = mean, verbose = false, cycle_lim = 30)
```

Cluster `input_vectors` using given distance and mean calculations, where a  vector with distance greater than `radius` from the nearest cluster forms a new cluster. Runs a maximum of `cycle_lim` iterations.

Returns (centroid_vectors, cluster_sizes, cluster_indices, centroid_indices) where cluster_indices is an array of arrays of indices of `input_vectors` grouped by cluster, and centroid_indices is an array of values where each value indexes the vector in `input_vectors` closest to the respective centroid.

This is a generalization of `dp_means`.


<a target='_blank' href='https://github.com/MurrellGroup/DPMeansClustering.jl/blob/25d48980e173119cfda767725fe058da77c60d99/src/dp.jl#L403-L416' class='documenter-source'>source</a><br>


```
dp_centers(inputs, radii::Vector{Float64}; distfunc = euclidean, center = mean, verbose = false, cycle_lim = [30...])
```

Cluster `input_vectors` using given distance and mean calculations recursively; `radii` should be an array of decreasing values, where each value is the radius of a successively finer clustering operation. During each clustering, a vector with distance greater than `radius` from the nearest cluster forms a new cluster. Runs a maximum of `cycle_lims` iterations during each cluster.

Returns (centroid_vectors, cluster_sizes, cluster_indices, centroid_indices) where cluster_indices is an array of arrays of indices of `input_vectors` grouped by cluster, and centroid_indices is an array of values where each value indexes the vector in `input_vectors` closest to the respective centroid.

This is a recursive implementation of `dp_centers(inputs, radius)` for faster clustering.


<a target='_blank' href='https://github.com/MurrellGroup/DPMeansClustering.jl/blob/25d48980e173119cfda767725fe058da77c60d99/src/dp.jl#L499-L513' class='documenter-source'>source</a><br>

<a id='DPMeansClustering.cluster' href='#DPMeansClustering.cluster'>#</a>
**`DPMeansClustering.cluster`** &mdash; *Function*.



```
cluster(inputs, radius::Float64; distfunc=euclidean, center=mean, verbose::Int64=1,
            cycle_lim::Int64=30, triangle=false)
```

A wrapper that decides whether or not to use triangle inequality.


<a target='_blank' href='https://github.com/MurrellGroup/DPMeansClustering.jl/blob/25d48980e173119cfda767725fe058da77c60d99/src/dp.jl#L2-L8' class='documenter-source'>source</a><br>

