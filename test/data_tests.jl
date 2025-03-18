longitudes =  [12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5]
latitudes = [-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5]

@testset "Test getMetaAttributesFromESMValToolConfigs" begin
end

@testset "Test addMetaData!" begin
end

@testset "Test getMetaDataFromYAML" begin
end

@testset "Test buildPathsToDataFiles" begin
end

@testset "Test getMetaDataID" begin
end

@testset "Test buildMetaData" begin
end

@testset "Test buildPathsForMetaAttrib" begin
end

@testset "Test applyDataConstraints!" begin
end

@testset "Test applyModelConstraints" begin
end

@testset "Test indexData" begin
end

@testset "Test getTimerangeAsAlias" begin
end

@testset "Test getAliasAsTimerange" begin
end

@testset "Test getRefPeriodAsTimerangeAndAlias" begin
end

@testset "Test subsetPaths" begin
end

@testset "Test subsetModelData" begin
end

@testset "Test loadDataFromMetadata" begin
end

@testset "Test getModelIDsFromPaths" begin
end

@testset "Test searchModelInPaths" begin
end

@testset "Test getSharedModelsFromPaths" begin
end

@testset "Test alignPhysics" begin
end

@testset "Test loadPreprocData" begin
end

@testset "Test loadDataFromESMValToolConfigs" begin
end

@testset "Test loadDataFromYAML" begin
end

@testset "Test averageEnsembleMembers!" begin
end

@testset "Test getPutDimArray" begin
    a = [1 2 3];
    b = [4 5 6];    
    da = DimArray(vcat(a,b), (Dim{:model}(["m1", "m2"]), Dim{:var}(["tas", "tos", "pr"])))
    @test mw.getAtModel(da, :model, "m1") == [1, 2, 3]

    mw.putAtModel!(da, :model, "m2", [17, 2, 1987])
    @test mw.getAtModel(da, :model, "m2") == [17, 2, 1987]
end



@testset "Test getGlobalMeans" begin
    a = [1.0 2.0 3.0 4.0];
    b = [4.0 5.0 6 5.0];
    da = YAXArray((Dim{:lon}(longitudes[1:2]), Dim{:lat}(latitudes[1:4])), vcat(a,b))
    gms = mw.getGlobalMeans(da)

    area_weights = mw.makeAreaWeightMatrix(Array(da.lon), Array(da.lat))
    aw_gm = sum(da .* area_weights)

end