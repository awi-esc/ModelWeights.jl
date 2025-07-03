import ModelWeights as mw
import ModelWeights.Data as mwd

using DimensionalData
using YAXArrays


longitudes =  [12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5]
latitudes = [-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5]
arr1 = YAXArray((Dim{:lat}(latitudes),Dim{:lon}(longitudes), Dim{:member}(["ESM1#r1i1p1f1", "ESM2#r1i1p1f1", "ESM2#r1i1p1f2"])), rand(7, 9, 3))
arr2 = YAXArray((Dim{:lat}(latitudes),Dim{:lon}(longitudes), Dim{:member}(["ESM1#r1i1p1f1", "ESM2#r1i1p1f1"])), rand(7, 9, 2)) 

a = [1.0 2.0 missing 4];
b = [5.0 missing missing 8];
c = [missing -1 0 1];

d1 = reshape(vcat([1 2], [3 4], [5 6]), 2, 3, 1);
d2 = reshape(1:6, 2, 3, 1);
yax1 = YAXArray((Dim{:row}(["r1", "r2"]),Dim{:column}(["c1", "c2", "c3"]), Dim{:stack}(["s1"])), d1)
yax2 = YAXArray((Dim{:row}(["r1", "r2"]),Dim{:column}(["c1", "c2", "c3"]), Dim{:stack}(["s1"])), d2) 


@testset "Test joinDataMaps" begin
    dm1 = mwd.DataMap(Dict("id1" => yax1))
    dm2 = mwd.DataMap(Dict("id2" => yax2))
    @test mwd.joinDataMaps(dm1, dm2) == mwd.DataMap(Dict(
        "id1" => yax1,
        "id2" => yax2,
    ))
    # warning is thrown and value of shared key is taken from the second argument
    @test (@test_logs (:warn,) mwd.joinDataMaps(dm1, mwd.DataMap(Dict("id1" => yax2)); warn_msg="my warning")) == mwd.DataMap(Dict("id1" => yax2))
end

@testset "Test metaDataFromESMValToolRecipes" begin
end

@testset "Test metaDataFromYAML" begin
end

@testset "Test resolvePathsFromMetaData" begin
end

@testset "Test filterPathsSharedModels" begin
    m1_tos = "lgm-cmip5-tos-climatologies/CMIP5_CNRM-CM5_Omon_lgm_r1i1p1_tos.nc"
    m1_tas = "lgm-cmip5-tas-climatologies/CMIP5_CNRM-CM5_Omon_lgm_r1i1p1_tas.nc"
    m2_tos = "lgm-cmip5-tos-climatologies/CMIP5_FGOALS-g2_Omon_lgm_r1i1p1_tos.nc"
    m2_tas = "lgm-cmip5-tas-climatologies/CMIP5_FGOALS-g2_Omon_lgm_r1i1p1_tas.nc"  
    
    m31_tas = "lgm-cmip5-tas-climatologies/CMIP5_GISS-E2-R_Amon_lgm_r1i1p150_tas.nc"
    m32_tas = "lgm-cmip5-tas-climatologies/CMIP5_GISS-E2-R_Amon_lgm_r1i1p151_tas.nc"
    
    m31_tos = "lgm-cmip5-tos-climatologies/CMIP5_GISS-E2-R_Amon_lgm_r1i1p150_tos.nc"
    m32_tos = "lgm-cmip5-tos-climatologies/CMIP5_GISS-E2-R_Amon_lgm_r1i1p151_tos.nc"
    
    paths_tas = [m1_tas, m2_tas, m31_tas, m32_tas]
    paths_tos = [m1_tos, m2_tos, m31_tos]
    paths_all = [paths_tos, paths_tas]
    n = length(vcat(paths_all...))
    fn_format = :esmvaltool
    paths_filtered = mwd.filterPathsSharedModels(paths_all, mwd.MODEL_LEVEL, fn_format)
    @test length(vcat(paths_filtered...)) == n

    paths_filtered = mwd.filterPathsSharedModels(paths_all, mwd.MEMBER_LEVEL, fn_format)
    @test !(m32_tas in paths_filtered[1])
    @test length(vcat(paths_filtered...)) == n-1
end


@testset "Test constrainMetaData!" begin
end


@testset "Test subsetModelData" begin
end


@testset "Test sharedModelsFromPaths" begin
end

@testset "Test alignPhysics" begin
end

@testset "Test loadPreprocData" begin
end

@testset "Test loadDataFromESMValToolRecipes" begin
end

@testset "Test loadDataFromYAML" begin
    # path = "/albedo/home/brgrus001/ModelWeightsPaper/work/configs/climwip/climwip-model-data.yml"
    # data = mwd.loadDataFromYAML(path, preview = true);
end

@testset "Test averageEnsembleMembers" begin
end

@testset "Test getAtModel" begin
    a = [1 2 3];
    b = [4 5 6];    
    da = DimArray(vcat(a,b), (Dim{:model}(["m1", "m2"]), Dim{:var}(["tas", "tos", "pr"])))
    @test mwd.getAtModel(da, :model, "m1") == [1, 2, 3]

    mwd.putAtModel!(da, :model, "m2", [17, 2, 1987])
    @test mwd.getAtModel(da, :model, "m2") == [17, 2, 1987]
end



@testset "Test defineDataMap" begin
    paths_lgm_tas = ["data/lgm-cmip5-tas-climatologies", "data/lgm-cmip6-tas-climatologies"]
    paths_lgm_tos = ["data/lgm-cmip5-tos-climatologies", "data/lgm-cmip6-tos-climatologies"]
    paths_lgm_tas = joinpath.(pwd(), paths_lgm_tas)
    paths_lgm_tos = joinpath.(pwd(), paths_lgm_tos)
    filename_format = :esmvaltool
    dtype = "cmip"

    paths_lgm = [paths_lgm_tas, paths_lgm_tos]
    ids = ["lgm-tas-data", "tos_lgm"]
    data = mw.defineDataMap(paths_lgm, ids; dtype, filename_format)
    @test size(data["lgm-tas-data"]) == (72, 36, 17)
    @test size(data["tos_lgm"]) == (72, 36, 16)

    constraint = Dict{String, Union{Vector{String}, Symbol}}("level_shared" => :member)
    data = mw.defineDataMap(paths_lgm, ids; dtype, filename_format, constraint)
    @test data["lgm-tas-data"].member == data["tos_lgm"].member
    @test length(data["lgm-tas-data"].member) == 15

    constraint["models"] = ["GISS-E2-R"]
    data = mw.defineDataMap(paths_lgm, ids; dtype, filename_format, constraint)
    @test length(data["lgm-tas-data"].member) == 2

    constraint["models"] = ["GISS-E2-R#r1i1p150"]
    data = mw.defineDataMap(paths_lgm, ids; dtype, filename_format, constraint)
    @test length(data["lgm-tas-data"].member) == 1
end

@testset "Test sharedLevels" begin
    data = mw.defineDataMap([arr1, arr2],  ["lgm-tos", "lgm-tas"])
    shared_members = mwd.sharedLevelMembers(data)
    @test all(x -> x in shared_members, ["ESM1#r1i1p1f1", "ESM2#r1i1p1f1"]) && !("ESM2#r1i1p1f2" in shared_members)
    shared_models = mwd.sharedLevelModels(data)
    @test all(x -> x in shared_models, ["ESM1", "ESM2"]) && all(x -> !(x in shared_models), ["ESM1#r1i1p1f1", "ESM2#r1i1p1f1", "ESM2#r1i1p1f2"])
end