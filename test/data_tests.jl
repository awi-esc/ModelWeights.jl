using DimensionalData
using YAXArrays
import ModelWeights as mw
longitudes =  [12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5]
latitudes = [-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5]


a = [1.0 2.0 missing 4];
b = [5.0 missing missing 8];
c = [missing -1 0 1];

d1 = reshape(vcat([1 2], [3 4], [5 6]), 2, 3, 1);
d2 = reshape(1:6, 2, 3, 1);
yax1 = YAXArray((Dim{:row}(["r1", "r2"]),Dim{:column}(["c1", "c2", "c3"]), Dim{:stack}(["s1"])), d1)
yax2 = YAXArray((Dim{:row}(["r1", "r2"]),Dim{:column}(["c1", "c2", "c3"]), Dim{:stack}(["s1"])), d2) 


@testset "Test joinDataMaps" begin
    dm1 = mw.DataMap(Dict("id1" => yax1))
    dm2 = mw.DataMap(Dict("id2" => yax2))
    @test mw.joinDataMaps(dm1, dm2) == mw.DataMap(Dict(
        "id1" => yax1,
        "id2" => yax2,
    ))
    # warning is thrown and value of shared key is taken from the second argument
    @test (@test_logs (:warn,) mw.joinDataMaps(dm1, mw.DataMap(Dict("id1" => yax2)))) == mw.DataMap(Dict("id1" => yax2))
end

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

@testset "Test loadDataFromESMValToolRecipes" begin
end

@testset "Test loadData" begin
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



@testset "Test computeGlobalMeans" begin
    a = [1.0 2.0 3.0 4.0];
    b = [4.0 5.0 6 5.0];
    da = YAXArray((Dim{:lon}(longitudes[1:2]), Dim{:lat}(latitudes[1:4])), vcat(a,b))
    gms = mw.computeGlobalMeans(da)

    area_weights = mw.makeAreaWeightMatrix(Array(da.lon), Array(da.lat))
    aw_gm = sum(da .* area_weights)

end