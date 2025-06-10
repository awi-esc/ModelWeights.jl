@testset "Test combineAll" begin
    # unequal number of entries (first more)
    v1 = ["tos", "tas"]; v2 = ["CLIM"];
    expected1 = ["tas_CLIM", "tos_CLIM"]
    result1 = sort(ModelWeights.Data.combineAll(v1, v2))
    @test result1 == expected1

    # one/two input vectors empty -> should throw warning
    msg = "At least one input vector is empty -> empty vector returned!"
    @test_logs (:warn, msg) ModelWeights.Data.combineAll(Vector{String}(), Vector{String}())
    @test_logs (:warn, msg) ModelWeights.Data.combineAll(["a"], Vector{String}())

    # unequal number of arguments (second more)
    v1 = ["pr", "psl"]; v2 = ["CLIM", "std", "CLIM-mon"]
    result2 = sort(ModelWeights.Data.combineAll(v1, v2))
    expected2 = ["pr_CLIM", "pr_CLIM-mon", "pr_std", "psl_CLIM", "psl_CLIM-mon", "psl_std"]
    @test result2 == expected2
end

@testset "Test computeAreaWeights" begin
end

@testset "Test individuatePath" begin
end