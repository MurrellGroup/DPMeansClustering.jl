using Distances

@testset "dp" begin
    #DONE
    @testset "dp_means" begin
        inputs = [[5, 0, 0, 1, 1, 5],
                  [5, 1, 0, 0, 1, 5],
                  [5, 0, 1, 0, 1, 5],
                  [0, 4, 6, 2, 0, 0],
                  [0, 4, 6, 1, 1, 0],
                  [0, 4, 6, 1, 0, 1]]
        radius = 3.0
        μs, sizes, indices, centroids = dp_means(inputs, radius)
        @test μs == [[5, 1/3, 1/3, 1/3, 1, 5],
                     [0, 4, 6, 4/3, 1/3, 1/3]]
        @test sizes == [3, 3]
        @test indices == [[1, 2, 3], [4, 5, 6]]
        @test centroids == [1, 4]
    end
    
    #DONE
    @testset "C2MC_dp_centers" begin
        inputs = [[5, 0, 0, 1, 1, 5],
                  [5, 1, 0, 0, 1, 5],
                  [5, 0, 1, 0, 1, 5],
                  [0, 4, 6, 2, 0, 0],
                  [0, 4, 6, 1, 1, 0],
                  [0, 4, 6, 1, 0, 1]]
        radius = 3.0
        μs, sizes, indices, centroids = C2MC_dp_centers(inputs, radius)
        @test μs == [[5, 1/3, 1/3, 1/3, 1, 5],
                     [0, 4, 6, 4/3, 1/3, 1/3]]
        @test sizes == [3, 3]
        @test indices == [[1, 2, 3], [4, 5, 6]]
        @test centroids == [1, 4]
    end
    
    @testset "dp_centers" begin
        inputs = [[5, 0, 0, 1, 1, 5],
                  [5, 1, 0, 0, 1, 5],
                  [5, 0, 1, 0, 1, 5],
                  [0, 4, 6, 2, 0, 0],
                  [0, 4, 6, 1, 1, 0],
                  [0, 4, 6, 1, 0, 1]]
        radii = 3.0
        μs, sizes, indices, centroids = dp_centers(inputs, radii)
        @test μs == [[5, 1/3, 1/3, 1/3, 1, 5],
                     [0, 4, 6, 4/3, 1/3, 1/3]]
        @test sizes == [3, 3]
        @test indices == [[1, 2, 3], [4, 5, 6]]
        @test centroids == [1, 4]
    end
    
    #recursive version of dp_centers
    @testset "dp_centers" begin
        inputs = [[5, 0, 0, 1, 1, 5],
                  [5, 1, 0, 0, 1, 5],
                  [5, 0, 1, 0, 1, 5],
                  [0, 4, 6, 2, 0, 0],
                  [0, 4, 6, 1, 1, 0],
                  [0, 4, 6, 1, 0, 1]]
        radii = [5.0, 3.0]
        μs, sizes, indices, centroids = dp_centers(inputs, radii)
        @test μs == [[5, 1/3, 1/3, 1/3, 1, 5],
                     [0, 4, 6, 4/3, 1/3, 1/3]]
        @test sizes == [3, 3]
        @test indices == [[1, 2, 3], [4, 5, 6]]
        @test centroids == [1, 4]
    end

end
