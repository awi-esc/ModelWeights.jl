using Logging
using ModelWeights
using Test


@testset "ModelWeights tests" begin
    @testset "Models tests" begin
        include("models_tests.jl")
    end
    @testset "Data tests" begin
        include("data_tests.jl")
    end
    @testset "Utils tests" begin
        include("utils_tests.jl")
    end
    @testset "Metadata tests" begin
        include("metadata_tests.jl")
    end
    @testset "Weights tests" begin
        include("weights_tests.jl")
    end
end
