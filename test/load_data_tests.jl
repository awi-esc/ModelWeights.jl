using DimensionalData
include("data.jl")

@testset "Testset CMIP Data" begin
    modelData =  SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, "CLIM", ["CMIP"]);
    @test length(keys(modelData)) == 3;
    @test size(modelData["tas"]) == (20, 19, 7)
    dimensions = DimensionalData.dims(modelData["tas"], :model);
    @test  Array(dimensions) == ["ACCESS1-0", "ACCESS1.3", "BNU-ESM", "CCSM4", "CCSM4", "CCSM4", "CCSM4"] # new code
    @warn "ACCESS1.3 should be ACCESS1-3 according to filenames (and content of metadata.yml), but the actual metadata says ACCESS1.3 instead."
end


@testset "Testset Observational Data" begin
    modelData =  SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, "CLIM", ["ERA5"]);
    @test length(keys(modelData)) == 3;
    @test size(modelData["tas"]) == (20, 19, 1)
    
    dimensions = DimensionalData.dims(modelData["tas"], :model);
    @test length(dimensions) == 1
    @test  dimensions[1] == "ERA5"
end