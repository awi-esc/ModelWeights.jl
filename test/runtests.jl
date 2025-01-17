using ModelWeights
using Test

@testset "ModelWeights tests" begin
    @testset "Models tests" begin
        include("models_tests.jl")
    end
    # @testset "Metadata tests" begin
    #     include("metadata_tests.jl")
    # end
end