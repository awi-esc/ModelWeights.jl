using SimilarityWeights
using Test

@testset "SimilarityWeights tests" begin
    @testset "Load data" begin
        include("load_data_tests.jl")
    end

    @testset "distance matrices" begin
        include("rmse_tests.jl")
    end
end