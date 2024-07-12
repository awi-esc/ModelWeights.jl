using SimilarityWeights
using Test

@testset "SimilarityWeights tests" begin
    @testset "Load data" begin
        include("load_data_tests.jl")
    end
end