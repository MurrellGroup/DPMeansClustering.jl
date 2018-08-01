module DPMeansClustering 
    
	using Reexport

    @reexport using StatsBase
    @reexport using Distances
    
	include("dp.jl")

    export dp_centers, dp_means, cluster, get_initial_centroids, get_initial_metacentroids, get_possible_centroids, C2MC_dp_centers
end
