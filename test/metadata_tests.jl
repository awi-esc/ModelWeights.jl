@testset "Test updateMetadataFromMultipleFiles!" begin
end

@testset "Test joinMetadata" begin
end

@testset "Test updateGroupedDataMetadata" begin
end

@testset "Test fixModelNamesMetadata" begin
end

@testset "Test uniqueMemberID" begin
end

@testset "Test getCMIPModelsKey" begin
    models = ["m1", "m2"]
    none_present = Dict("random-key" => models)
    cmip5 = Dict("random_key"=>"random_value", "model_id" => models)
    cmip6 = Dict("source_id" => models)
    both = Dict("source_id" => models, "model_id" => models, "other_key" => ["other_val1", "other_val2"])

    @test ModelWeights.Data.getCMIPModelsKey(cmip5) == "model_id"
    @test ModelWeights.Data.getCMIPModelsKey(cmip6) == "source_id"
    
    @test ModelWeights.Data.getCMIPModelsKey(both) == "source_id"
    warning_both = "Dictionary contains keys source_id and model_id, source_id is used!"
    # TODO: check how this works!!
    # @test_logs (:debug, warning_both) ModelWeights.getCMIPModelsKey(both)
    
    @test_throws ArgumentError ModelWeights.Data.getCMIPModelsKey(none_present)
end