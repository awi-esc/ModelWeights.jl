@testset "Test updateMetadata!" begin
end

@testset "Test joinMetadata" begin
end

@testset "Test updateGroupedDataMetadata" begin
end

@testset "Test fixModelNamesMetadata" begin
end

@testset "Test getUniqueMemberIds" begin
end

@testset "Test getCMIPModelsKey" begin
    models = ["m1", "m2"]
    none_present = Dict("random-key" => models)
    cmip5 = Dict("random_key"=>"random_value", "model_id" => models)
    cmip6 = Dict("source_id" => models)
    both = Dict("source_id" => models, "model_id" => models, "other_key" => ["other_val1", "other_val2"])

    @test ModelWeights.getCMIPModelsKey(cmip5) == "model_id"
    @test ModelWeights.getCMIPModelsKey(cmip6) == "source_id"
    
    @test ModelWeights.getCMIPModelsKey(both) == "source_id"
    warning_both = "Dictionary contains  keys 'source_id' (CMIP6) and 'model_id' (CMIP5). 'source_id' is used!"
    @test_logs (:warn, warning_both) ModelWeights.getCMIPModelsKey(both)
    
    @test_throws ArgumentError ModelWeights.getCMIPModelsKey(none_present)
end