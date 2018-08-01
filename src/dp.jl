
"""
    cluster(inputs, radius::Float64; distfunc=euclidean, center=mean, verbose::Int64=1,
                cycle_lim::Int64=30, triangle=false)

A wrapper that decides whether or not to use triangle inequality.

"""
function cluster(inputs, radius::Float64; distfunc=euclidean, center=mean, verbose::Int64=1,
                cycle_lim::Int64=30, triangle=false)
    if triangle
        #use triangle, since euclid is a metric or approximation is okay
        return C2MC_dp_centers(inputs, radius; distfunc=distfunc, center=center, verbose=verbose,
                        cycle_lim=cycle_lim)
    else
        if verbose > 1
            verb = true
        else
            verb = false
        end
        #distfunc is non-euclidean and approximations are not acceptable
        return dp_centers(inputs, radius, distfunc=distfunc, center=center, verbose=verb,
                        cycle_lims=cycle_lim)
    end
end


"""
    get_initial_centroids(inputs, radius::Float64, input_length::Int64, distfunc, center, verbose::Int64)

Performs the first cycle of clustering in preparation for the creation of metacentroids
"""
function get_initial_centroids(inputs, radius::Float64, input_length::Int64, distfunc, center, verbose::Int64)
    num_C = 1
    C_arr = [inputs[1]]
    C_matrix = [falses(1, input_length)]
    for i in 1:length(inputs)
        curr_input = inputs[i]
        dists = [distfunc(C_arr[j], curr_input) for j in 1:num_C]
        if minimum(dists) > radius
            # Make a new cluster
            for k in 1:num_C #eliminate the input from its old cluster
                C_matrix[k][i] = false
            end
            num_C += 1 #increment clusters
            push!(C_matrix, falses(1, input_length)) #add to matrix
            C_matrix[num_C][i] = true #mark which input is in the new cluster
            push!(C_arr, curr_input)
        else
            # Assign point to best cluster
            min_index = indmin(dists)
            if !(C_matrix[min_index][i])
                for k in 1:num_C #eliminate the input from its old cluster
                    C_matrix[k][i] = false
                end
                C_matrix[min_index][i] = true #mark which input is in the new cluster
            end
        end
    end
    C_arr = [center(inputs[find(C_matrix[i])]) for i in 1:num_C]
    if verbose == 2
        println("cycle 1. n_clusters=$(num_C)")
    end
    return C_arr, C_matrix
end

"""
    get_initial_metacentroids(inputs, C_arr, MC_radius::Float64, distfunc, center, verbose::Int64)

Create metacentroids for the triangle inequality
"""
function get_initial_metacentroids(inputs, C_arr, MC_radius::Float64, distfunc, center, verbose::Int64)
    num_MC = 1
    changed = true
    num_C = length(C_arr)
    MC_arr = [C_arr[1]]
    MC_matrix = [falses(1, num_C)]
    C2MC_dist_matrix = Array{Array{Float64, 1}, 1}()
    changed = true
    for i in 1:num_C
        curr_C = C_arr[i]
        dists = [distfunc(MC_arr[j], curr_C) for j in 1:num_MC]
        if minimum(dists) > MC_radius
            for k in 1:num_MC #eliminate the input from its old cluster
                MC_matrix[k][i] = false
            end
            num_MC += 1 #increment clusters
            push!(MC_matrix, falses(1, num_C)) #add to matrix
            MC_matrix[num_MC][i] = true #mark which input is in the new cluster
            push!(MC_arr, curr_C)
        else
            # Assign point to best cluster
            min_index = indmin(dists)
            if !(MC_matrix[min_index][i])
                for k in 1:num_MC #eliminate the input from its old cluster
                    MC_matrix[k][i] = false
                end
                MC_matrix[min_index][i] = true #mark which input is in the new cluster
            end
        end
    end
    MC_arr = [center(C_arr[find(MC_matrix[MC_index])]) for MC_index in 1:num_MC]
    for MC_index in 1:num_MC #get distances between centroids and metacentroids
        C_indices = find(MC_matrix[MC_index])
        push!(C2MC_dist_matrix, [distfunc(MC_arr[MC_index], C_arr[index]) for index
                               in C_indices])
    end
    return MC_arr, C2MC_dist_matrix, MC_matrix
end

"""
    get_possible_centroids(MC_matrix, C_arr, R2MC_dist_arr::Array{Float64, 1},
                            C2MC_dist_matrix::Array{Array{Float64, 1}, 1}, radius::Float64, num_MC::Int64,
                            distfunc, MC_radius::Float64, verbose::Int64)

Uses the triangle inequality in order to determine which centroids a read may cluster to. This is
done by finding the difference between the read to metacentroid(R2MC_dist_arr) as well as centroid
to metacentroid distance(C2MC_dist_matrix), since the triangle inequality mandates that the
distance from read to centroid must be at least the difference between the other two sides

Returns (possible_centroids_indices) which contains the indices of the centroids a read may belong
to. The indices correspond to the centroids' locations in the C_arr.
"""

function get_possible_centroids(MC_matrix, C_arr, R2MC_dist_arr::Array{Float64, 1},
                            C2MC_dist_matrix::Array{Array{Float64, 1}, 1}, radius::Float64, num_MC::Int64,
                            distfunc, MC_radius::Float64, verbose::Int64)
    possible_centroids_indices = Array{Int32, 1}()
    total_radius = MC_radius + radius
    possible_metacentroids_indices = Array{Int32, 1}()
    for MC_index in 1:num_MC
        R2MC = R2MC_dist_arr[MC_index]
        if R2MC <= total_radius  #only reads at a certain dist from an MC qualifies
            push!(possible_metacentroids_indices, MC_index)
        end
    end

    for MC in possible_metacentroids_indices #go through each predetermined MC
        R2MC_arr = R2MC_dist_arr[MC]
        C_indices = find(MC_matrix[MC]) #indices of centroids belonging to this MC
        C2MC_distances = C2MC_dist_matrix[MC] #distances between each centroid of this MC
        num_possible_C = length(C2MC_distances)
        arr = abs.(C2MC_distances - R2MC_arr) #minimum read to centroid distance
        for C_index in 1:num_possible_C # go through each C in this MC
            if(arr[C_index] <= radius) #the minimum read to centroid must be less than radius
                push!(possible_centroids_indices, C_indices[C_index])
            end
        end
    end
    return possible_centroids_indices
end

"""
    C2MC_dp_centers(inputs, radius::Float64; distfunc=euclidean, center=mean, verbose::Int64=1,
                        cycle_lim::Int64=30)

Clusters the input kmer arrays using the triangle inequality for given radius
"""
function C2MC_dp_centers(inputs, radius::Float64; distfunc=euclidean, center=mean, verbose::Int64=1,
                        cycle_lim::Int64=30)
    input_length = length(inputs)

    #1. Get initial centroids. Looks correct(ran comparison with old_dp)
    C_arr, C_matrix = get_initial_centroids(inputs, radius, input_length, distfunc,
                                            center, verbose)
    MC_radius =  .5 * radius #how to determine optimum radius size for MCs?
    num_C = length(C_arr)

    #2. Create initial metacentroids as well as C2MC distances
    MC_arr, C2MC_dist_matrix, MC_matrix = get_initial_metacentroids(inputs,C_arr, MC_radius, distfunc,
                                                                    center, verbose)
    num_MC = length(MC_arr)
    #3. Get a matrix of the distances between the reads and metacentroids
    R2MC_dist_matrix = [distfunc(MC_arr[MC_index], inputs[input_index]) for input_index in
                        1:input_length, MC_index in 1:num_MC]

    total_proportions = 0
    cycles = 1 #initial centroid cycle already complete
    changed = true

    while changed && cycles < cycle_lim
        changed = false
        for i in 1:input_length
            curr_input = inputs[i]
            R2MC_dist_arr = R2MC_dist_matrix[i,:]
            #5. Triangle to get possible C's
            possible_centroids_arr = get_possible_centroids(MC_matrix, C_arr, R2MC_dist_arr,
                                                C2MC_dist_matrix, radius, num_MC, distfunc, MC_radius, verbose)
            dists = [distfunc(C_arr[j], curr_input) for j in possible_centroids_arr]

            new_proportion = length(possible_centroids_arr) / num_C
            total_proportions+=new_proportion
            #6. Insert read into C
            #the read neads a new centroid
            if length(dists) == 0 || minimum(dists) > radius
                # Make a new cluster
                for k in 1:num_C #eliminate the input from its old cluster
                    C_matrix[k][i] = false
                end
                num_C += 1 #increment clusters
                push!(C_matrix, falses(1, input_length)) #add to matrix
                C_matrix[num_C][i] = true #mark which input is in the new cluster
                push!(C_arr, curr_input)
                changed = true

                MC_distances = [distfunc(MC_arr[MC_index], curr_input) for MC_index in 1:num_MC]
                min_MC_index = indmin(MC_distances)
                minimum_distance = minimum(MC_distances)
                #the new centroid can be added to an existing metacentroid
                if minimum_distance <= MC_radius + radius
                    #add C to MC_matrix
                    for i in 1:num_MC
                        MC_matrix[i] = hcat(MC_matrix[i], false)
                    end
                    MC_matrix[min_MC_index][num_C] = true
                    push!(C2MC_dist_matrix[min_MC_index], minimum_distance)
                else #the new centroid requires a new metacentroid
                    for i in 1:num_MC
                        MC_matrix[i] = hcat(MC_matrix[i], false)
                    end
                    num_MC+=1

                    push!(MC_matrix, falses(1, num_C))
                    MC_matrix[num_MC][num_C] = true
                    push!(MC_arr, curr_input)
                    push!(C2MC_dist_matrix, [0.0])
                    new_dist_arr = [distfunc(inputs[j], curr_input) for j in 1:input_length]
                    R2MC_dist_matrix = hcat(R2MC_dist_matrix, new_dist_arr)
                end
            else
                # Assign point to best cluster
                min_index = possible_centroids_arr[indmin(dists)]
                if !(C_matrix[min_index][i]) #do a check to see if already in this cluster
                    for k in 1:num_C #eliminate the input from its old cluster
                        C_matrix[k][i] = false
                    end
                    C_matrix[min_index][i] = true #mark which cluster the input was moved to
                    changed= true
                end
            end
        end

        #7. Update C's
        C_arr = [center(inputs[find(C_matrix[i])]) for i in 1:num_C]

        #8. Update C2MC. Exactly the same as getting initials.
        cycles+=1

        empty!(C2MC_dist_matrix)
        for MC_index in 1:num_MC #get distances between centroids and metacentroids
            C_indices = find(MC_matrix[MC_index])
            push!(C2MC_dist_matrix, [distfunc(MC_arr[MC_index], C_arr[index]) for index
                                   in C_indices])
        end
        if verbose == 2
            println("cycle $(cycles). n_clusters=$(num_C)")
        end
    end

    total_proportions = total_proportions / (input_length * (cycles - 1))
    if verbose > 0
        println("A clustering call has converged!")
        println("The average proportion of centroids passing triangle inequality
                is $(total_proportions)")
        println("There are $(num_C) centroids and $(num_MC) metacentroids")
    end

    cluster_indices = [find(C_matrix[i]) for i in 1:num_C]
    # Cleaning out the empty clusters. Maybe only run this if there
    # are length 0 elements in cluster_indices.
    sizes = zeros(Int, num_C)
    for i in 1:num_C
        sizes[i] = length(cluster_indices[i])
    end

    clean_cluster_indices = Vector{Int}[]
    clean_μs = []  # TODO: set  type. But mean may not be same type as inputs!!!
    clean_sizes = Int[]
    for i in 1:num_C
        if length(cluster_indices[i]) > 0
            push!(clean_cluster_indices, cluster_indices[i])
            push!(clean_μs, C_arr[i])
            push!(clean_sizes, sizes[i])
        end
    end

    n_keep = length(clean_μs)
    centroid_inds = [-1 for i in 1:n_keep]
    for j in 1:n_keep
    dists = [distfunc(clean_μs[j], inputs[clean_cluster_indices[j][i]])
    for i in 1:length(clean_cluster_indices[j])]
        centroid_inds[j] = clean_cluster_indices[j][indmin(dists)]
    end
    return clean_μs, clean_sizes, clean_cluster_indices, centroid_inds
end

"""
    dp_means(input_vectors::Array{Vector{Int}, 1}, radius::Float64; verbose = false)

Cluster `input_vectors` using euclidean distance metric and arithmetic mean, where a 
vector with distance greater than `radius` from the nearest cluster forms a new cluster.

Returns (centroid_vectors, cluster_sizes, cluster_indices, centroid_indices) where cluster_indices
is an array of arrays of indices of `input_vectors` grouped by cluster, and centroid_indices
is an array of values where each value indexes the vector in `input_vectors` closest to
the respective centroid.

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
"""
function dp_means(input_vectors::Vector{Vector{Int}}, radius::Float64; verbose = false)
    # TODO:
    # 0) Return exemplar indices for each cluster, where the vector
    #    closest to the cluster mean is selected. [DONE]
    # 0.1) Check if exemplars close to cluster means tend to have better accuracies.
    # 0.2) Return cluster assignments more conveniently. [DONE]
    # 1) Experiment with removing singleton clusters after convergence
    #    and starting the iteration again. [maybe nah]
    n_clusters = 1
    oldZs = [-1 for i in 1:length(input_vectors)]
    μs = [input_vectors[1]]
    newZs = [1 for i in 1:length(input_vectors)]
    sizes = zeros(n_clusters)
    cycles = 0
    while oldZs != newZs
        oldZs = copy(newZs)
        for i in 1:length(input_vectors)
            dists = [euclidean(μs[j], input_vectors[i]) for j in 1:n_clusters]
            if minimum(dists) > radius
                # Make a new cluster
                n_clusters += 1
                newZs[i] = n_clusters  #  Note we've just incremented n_clusters!
                push!(μs, input_vectors[i])
            else
                # Assign point to best cluster
                newZs[i] = indmin(dists)
            end
        end
        # Reset all means
        μs = μs * 0.0
        sizes = zeros(n_clusters)
        for i in 1:length(input_vectors)
            μs[newZs[i]] += input_vectors[i]
            sizes[newZs[i]] += 1
        end
        for j in 1:n_clusters
            μs[j] = μs[j] / sizes[j]
        end
        cycles += 1
        if verbose
            println("cycle $(cycles). n_clusters=$(n_clusters)")
        end
    end
    if verbose
        println("Converged!")
    end
    cluster_indices = [[] for i in 1:n_clusters]
    for i in 1:length(input_vectors)
        push!(cluster_indices[newZs[i]], i)
    end

    # Cleaning out the empty clusters. Maybe only run this if there
    # are length 0 elements in cluster_indices.
    clean_cluster_indices = Vector{Int}[]
    clean_μs = Vector{Float64}[]
    clean_sizes = Int[]
    for i in 1:length(cluster_indices)
        if length(cluster_indices[i]) > 0
            push!(clean_cluster_indices, cluster_indices[i])
            push!(clean_μs, μs[i])
            push!(clean_sizes, sizes[i])
        end
    end

    n_keep = length(clean_μs)
    if verbose
        println("keeping $(n_keep) after removing empty clusters.")
    end
    centroid_indices = [-1 for i in 1:n_keep]
    for j in 1:n_keep
        dists = [euclidean(clean_μs[j], input_vectors[clean_cluster_indices[j][i]])
                 for i in 1:length(clean_cluster_indices[j])]
        centroid_indices[j] = clean_cluster_indices[j][indmin(dists)]
    end
    return clean_μs, clean_sizes, clean_cluster_indices, centroid_indices
end

"""
    dp_centers(inputs, radius::Float64; distfunc = euclidean, center = mean, verbose = false, cycle_lim = 30)

Cluster `input_vectors` using given distance and mean calculations, where a 
vector with distance greater than `radius` from the nearest cluster forms a new cluster.
Runs a maximum of `cycle_lim` iterations.

Returns (centroid_vectors, cluster_sizes, cluster_indices, centroid_indices) where cluster_indices
is an array of arrays of indices of `input_vectors` grouped by cluster, and centroid_indices
is an array of values where each value indexes the vector in `input_vectors` closest to
the respective centroid.

This is a generalization of `dp_means`.
"""
function dp_centers(inputs, radius::Float64;
                    distfunc = euclidean,
                    center = mean,
                    verbose = false, 
                    cycle_lims = 30)

    # ToDo:
    # 0) Return exemplar indices for each cluster, where the vector
    #    closest to the cluster mean is selected. [DONE]
    # 0.1) Check if exemplars close to cluster means tend to have better accuracies. [YES THEY DOOOOOOO!]
    # 0.2) Return cluster assignments more conveniently. [DONE]
    # 1) Experiment with removing singleton clusters after convergence
    #    and starting the iteration again.
    n_clusters = 1
    oldZs = [-1 for i in 1:length(inputs)]
    μs = [inputs[1]]
    newZs = [1 for i in 1:length(inputs)]
    sizes = zeros(Int, n_clusters)
    cluster_indices = [[] for i in 1:n_clusters]
    cycles = 0
    while oldZs != newZs && cycles < cycle_lims
        oldZs = copy(newZs)
        for i in 1:length(inputs)
            dists = [distfunc(μs[j], inputs[i]) for j in 1:n_clusters]
            if minimum(dists) > radius
                # Make a new cluster
                n_clusters += 1
                # Note we've just incremented n_clusters!
                newZs[i] = n_clusters
                push!(μs, inputs[i])
            else
                # Assign point to best cluster
                newZs[i] = indmin(dists)
            end
        end

        # Reset all means
        cluster_indices = [[] for i in 1:n_clusters]
        for i in 1:length(inputs)
            push!(cluster_indices[newZs[i]], i)
        end

        μs = [center(inputs[cluster_indices[i]]) for i in 1:n_clusters]
        cycles += 1
        if verbose
            println("cycle $(cycles). n_clusters=$(n_clusters)")
        end
    end
    if verbose
        println("Converged!")
    end
    # Cleaning out the empty clusters. Maybe only run this if there
    # are length 0 elements in cluster_indices.
    sizes = zeros(Int, n_clusters)
    for i in 1:length(inputs)
        sizes[newZs[i]] += 1
    end
    clean_cluster_indices = Vector{Int}[]
    clean_μs = []  # TODO: set  type. But mean may not be same type as inputs!!!
    clean_sizes = Int[]
    for i in 1:length(cluster_indices)
        if length(cluster_indices[i]) > 0
            push!(clean_cluster_indices, cluster_indices[i])
            push!(clean_μs, μs[i])
            push!(clean_sizes, sizes[i])
        end
    end

    n_keep = length(clean_μs)
    if verbose
        println("n_keep=$(n_keep) after removing empty clusters.")
    end

    centroid_inds = [-1 for i in 1:n_keep]
    for j in 1:n_keep
        dists = [distfunc(clean_μs[j], inputs[clean_cluster_indices[j][i]])
                 for i in 1:length(clean_cluster_indices[j])]
        centroid_inds[j] = clean_cluster_indices[j][indmin(dists)]
    end
    return clean_μs, clean_sizes, clean_cluster_indices, centroid_inds
end

"""
    dp_centers(inputs, radii::Vector{Float64}; distfunc = euclidean, center = mean, verbose = false, cycle_lim = [30...])

Cluster `input_vectors` using given distance and mean calculations recursively; `radii` should
be an array of decreasing values, where each value is the radius of a successively finer clustering
operation. During each clustering, a vector with distance greater than `radius` from the nearest cluster forms a new cluster.
Runs a maximum of `cycle_lims` iterations during each cluster.

Returns (centroid_vectors, cluster_sizes, cluster_indices, centroid_indices) where cluster_indices
is an array of arrays of indices of `input_vectors` grouped by cluster, and centroid_indices
is an array of values where each value indexes the vector in `input_vectors` closest to
the respective centroid.

This is a recursive implementation of `dp_centers(inputs, radius)` for faster clustering.
"""
function dp_centers(inputs, radii::Array{Float64};
                    distfunc = euclidean,
                    center = mean,
                    verbose = false, 
                    cycle_lims = 20)#fill(30, length(radii)))
    cycle_lims = 20; #remove this later
    μs, sizes, cluster_inds, centroid_inds = dp_centers(inputs, radii[1]; 
                                                        distfunc=distfunc,
                                                        center=center,
                                                        verbose=verbose, 
                                                        cycle_lims=cycle_lims)#cycle_lims[1])
    if length(radii) > 1
        newClustsZ = [collect(dp_centers(inputs[cluster_inds[i]], radii[2:end],
                                                         distfunc=distfunc, 
                                                         center=center, 
                                                         verbose=verbose, 
                                                         cycle_lims=cycle_lims))#cycle_lims[2:end]))
                for i in 1:length(cluster_inds)]

        # reformat to be consistent with nonrecursive dp_centers output.
        # "transpose" arrays, ie. [[[a, a, a], [b, b, b]], [[A, A, A], [B, B, B]]] to
        # [[[a, a, a], [A, A, A]], [[b, b, b], [B, B, B]]]
        newClusts = [[newClustsZ[i][j] for i in 1:length(newClustsZ)] for j in 1:4]

        # new, finer sizes, means, and centroid inds are calculated
        # but here we get the new cluster_inds and centroid_inds in terms of the newly computed clusters
        for i in 1:length(cluster_inds)
            for j in 1:length(newClusts[3][i])
                newClusts[3][i][j] = cluster_inds[i][newClusts[3][i][j]]
            end
            # merge finer clusters into bigger one to get relative indices of new centroids
            # FIX: currently relies on dp_centers returning sorted lists of indices
            newClusts[4][i] = sort(vcat(newClusts[3][i]...))[newClusts[4][i]]
        end
        
        # merge sub-clusters into a single array
        for i in 1:4
            newClusts[i] = vcat(newClusts[i]...)
        end
        return newClusts
    else
        return μs, sizes, cluster_inds, centroid_inds
    end
end


