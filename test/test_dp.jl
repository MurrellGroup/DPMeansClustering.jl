using Distances

@testset "dp" begin
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

    @testset "dp_centers" begin
        inputs = [[5, 0, 0, 1, 1, 5],
                  [5, 1, 0, 0, 1, 5],
                  [5, 0, 1, 0, 1, 5],
                  [0, 4, 6, 2, 0, 0],
                  [0, 4, 6, 1, 1, 0],
                  [0, 4, 6, 1, 0, 1]]
        radius = 3.0
        μs, sizes, indices, centroids = dp_centers(inputs, radius)
        @test μs == [[5, 1/3, 1/3, 1/3, 1, 5],
                     [0, 4, 6, 4/3, 1/3, 1/3]]
        @test sizes == [3, 3]
        @test indices == [[1, 2, 3], [4, 5, 6]]
        @test centroids == [1, 4]
    end

    @testset "dp_centers with tuples" begin
        inputs = [([5, 0], "ACGT"),
                  ([4, 1], "ACGT"),
                  ([5, 0], "ACGG"),
                  ([0, 5], "TGCA"),
                  ([1, 4], "TGCA"),
                  ([0, 5], "TGCC")]
        radius = 3.0

        function distfunc(x, y)
            return euclidean(x[1], y[1]) + sum(collect(x[2]) .!= collect(y[2]))
        end

        function center(xs)
            chars = hcat([collect(x[2]) for x in xs]...)
            return (mean([x[1] for x in xs]), join([mode(chars[i,:]) for i in 1:size(chars)[1]]))
        end

        μs, sizes, indices, centroids = dp_centers(inputs, radius;
                                                   distfunc=distfunc, center=center)

        @test μs == [([14/3, 1/3], "ACGT"),
                     ([1/3, 14/3], "TGCA")]
        @test sizes == [3, 3]
        @test indices == [[1, 2, 3], [4, 5, 6]]
        @test centroids == [1, 4]
    end
end
