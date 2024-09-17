using DimensionalData
include("data.jl")

@testset "Testset CMIP Data" begin
    modelData =  SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["CMIP"]);
    @test length(keys(modelData["CLIM"])) == 3;
    @test size(modelData["CLIM"]["tas"]) == (20, 19, 7)
    dimensions = DimensionalData.dims(modelData["CLIM"]["tas"], :model);
    @test  Array(dimensions) == ["ACCESS1-0", "ACCESS1.3", "BNU-ESM", "CCSM4", "CCSM4", "CCSM4", "CCSM4"]
    @warn "ACCESS1.3 should be ACCESS1-3 according to filenames (and content of metadata.yml), but the actual metadata says ACCESS1.3 instead."
end


@testset "Testset Observational Data" begin
    modelData =  SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["ERA5"]);
    @test length(keys(modelData["CLIM"])) == 3;
    @test size(modelData["CLIM"]["tas"]) == (20, 19, 1)
    
    dimensions = DimensionalData.dims(modelData["CLIM"]["tas"], :model);
    @test length(dimensions) == 1
    @test  dimensions[1] == "native6_ERA5_reanaly_v1_Amon_tas_1995-2014"
end