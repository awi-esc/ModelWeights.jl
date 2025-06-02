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
    # path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical"
    # paths_cems2 = mw.buildPathsToDataFiles(path_data, true; model_constraints=["CESM2"])
end

@testset "Test getMetaDataID" begin
end

@testset "Test buildMetaData" begin
end

@testset "Test buildPathsForMetaAttrib" begin
    base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical"
    attribs = Dict(
        "_variable" => "tas",
        "_statistic" => "CLIM",
        "_alias" => "historical",
    )
    dir_per_var = true
    paths = mw.buildPathsForMetaAttrib(base_path, attribs, dir_per_var)
    subdir = "recipe_cmip5_historical_tas_20250211_094633"
    diagnostic = attribs["_variable"] * "_" * attribs["_statistic"]
    p1 = joinpath(base_path, subdir, "preproc", attribs["_alias"], diagnostic)
    @assert isdir(p1)
    @test p1 in paths

    p2 = joinpath(base_path, subdir, "preprocessor", attribs["_alias"], diagnostic)
    @assert !isdir(p2)
    @test !(p2 in paths)

    subdir = "recipe_cmip5_historical_psl_timeseries_20250228_084113"
    p3 = joinpath(base_path, subdir, "preproc", attribs["_alias"], "psl_CLIM-ann")
    @assert isdir(p3)
    @test !(p3 in paths)

    subdir = "recipe_cmip6_historical_tas_20250207_080843"
    p4 = joinpath(base_path, subdir, "preproc", attribs["_alias"], diagnostic)
    @assert isdir(p4)
    @test p4 in paths
    paths = mw.buildPathsForMetaAttrib(
        base_path, attribs, dir_per_var; 
        subdir_constraints = ["20250211"]
    )
    @test !(p4 in paths)
    @test p1 in paths
end


@testset "Test getPathsToData" begin
    path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical"
    meta = Dict(
        "_variable" => "tas",
        "_statistic" => "CLIM",
        "_alias" => "historical"
    )
    dir_per_var = true
    is_model_data = true
    ds_constraint = Dict(
        "models" =>  ["ACCESS-CM2", "CESM2#r1i1p1f1"]
    )
    paths_to_files = mw.getPathsToData(
        meta, path_data, dir_per_var, is_model_data; constraint = ds_constraint,
    )
    base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical"
    subdir = "recipe_cmip6_historical_tas_20250207_080843"
    alias = "historical"
    diagnostic = meta["_variable"] * "_" * meta["_statistic"]
    p_access = joinpath(base_path, subdir, "preproc", alias, diagnostic, "CMIP6_ACCESS-CM2_Amon_historical_r1i1p1f1_tas_gn.nc")
    p_cesm2 = joinpath(base_path, subdir, "preproc", alias, diagnostic, "CMIP6_CESM2_Amon_historical_r1i1p1f1_tas_gn.nc")
    p_other = joinpath(base_path, subdir, "preproc", alias, diagnostic, "CMIP6_BCC-CSM2-MR_Amon_historical_r2i1p1f1_tas_gn.nc")
    @assert isfile(p_access) && isfile(p_cesm2) && isfile(p_other)
    @test p_access in paths_to_files
    @test p_cesm2 in paths_to_files
    @test !(p_other in paths_to_files)
end


@testset "Test filterPathsSharedModels!" begin
    p1_tos = "LGM/recipe_cmip5_lgm_tos_20241114_150049/preproc/lgm/tos_CLIM/CMIP5_CNRM-CM5_Omon_lgm_r1i1p1_tos.nc"
    p2_tos = "LGM/recipe_cmip5_lgm_tos_20241114_150049/preproc/lgm/tos_CLIM/CMIP5_FGOALS-g2_Omon_lgm_r1i1p1_tos.nc"
    p1_tas = "LGM/recipe_cmip5_lgm_tas_20241114_150049/preproc/lgm/tas_CLIM/CMIP5_CNRM-CM5_Omon_lgm_r1i1p1_tas.nc"
    p2_tas = "LGM/recipe_cmip5_lgm_tas_20241114_150049/preproc/lgm/tas_CLIM/CMIP5_FGOALS-g2_Omon_lgm_r1i1p1_tas.nc"
    meta = Dict(
        "tos_CLIM_historical" => Dict("_paths" => [p1_tos, p2_tos]),   
        "tas_CLIM_historical" => Dict("_paths" => [replace(p1_tas, "r1i1p1"=>"r2i1p1"), p2_tas])     
    )
    mw.filterPathsSharedModels!(meta, mw.MODEL)
    @test all(x -> x in meta["tos_CLIM_historical"]["_paths"], [p1_tos, p2_tos])
    mw.filterPathsSharedModels!(meta, mw.MEMBER)
    @test p2_tos in meta["tos_CLIM_historical"]["_paths"]
    @test !(p1_tos  in meta["tos_CLIM_historical"]["_paths"])
end


@testset "Test applyDataConstraints!" begin
end

@testset "Test applyModelConstraints" begin
    path_no = "path/to/ACCESS_r1i1p1f1_gr.nc" # does not stick to name convention (model must be in between two underscores (_model_)
    path_ok = "cmip_ACCESS_r1i1p1f1_gr_1970.nc"
    path_ok_notime = "cmip_ACCESS_r1i1p1f1_gr.nc"
    path_other_member = "cmip_ACCESS_r1i1p1f2_gr.nc"

    allowed_models = ["ACCESS#r1i1p1f1"]
    @test !mw.applyModelConstraints(path_no, allowed_models)
    @test mw.applyModelConstraints(path_ok, allowed_models)
    @test mw.applyModelConstraints(path_ok_notime, allowed_models)
    @test !mw.applyModelConstraints(path_other_member, allowed_models)

    allowed_models = ["ACCESS"]
    @test mw.applyModelConstraints(path_other_member, allowed_models)

    allowed_models = ["ACCESS#r1i1p1f1_gr"]
    @test !mw.applyModelConstraints("cmip_ACCESS_r1i1p1f1_gn.nc", allowed_models)
    @test mw.applyModelConstraints("cmip_ACCESS_r1i1p1f1_gr.nc", allowed_models)
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
    member_present = "CESM2#r1i1p1f1"
    model_present = "CESM2"
    model_absent = "AWI-ESM"
    paths = [
        "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_cmip6_historical_tas_20250207_080843/preproc/historical/tas_CLIM/CMIP6_ACCESS-CM2_Amon_historical_r1i1p1f1_tas_gn.nc",
        "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_cmip6_historical_tas_20250207_080843/preproc/historical/tas_CLIM/CMIP6_CESM2_Amon_historical_r1i1p1f1_tas_gn.nc"
    ]
    @test mw.searchModelInPaths(member_present, paths)
    @test mw.searchModelInPaths(model_present, paths)
    @test !mw.searchModelInPaths(model_absent, paths)
end

@testset "Test getSharedModelsFromPaths" begin
end

@testset "Test alignPhysics" begin
end

@testset "Test loadPreprocData" begin
end

@testset "Test loadDataFromESMValToolRecipes" begin
end

@testset "Test loadDataFromYAML" begin
    # path = "/albedo/home/brgrus001/ModelWeightsPaper/work/configs/climwip/climwip-model-data.yml"
    # data = mw.loadDataFromYAML(path, preview = true);
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



# @testset "Test computeGlobalMeans" begin
#     a = [1.0 2.0 3.0 4.0];
#     b = [4.0 5.0 6 5.0];
#     da = YAXArray((Dim{:lon}(longitudes[1:2]), Dim{:lat}(latitudes[1:4])), vcat(a,b))
#     gms = mw.computeGlobalMeans(da)

#     area_weights = mw.makeAreaWeightMatrix(Array(da.lon), Array(da.lat))
#     aw_gm = sum(da .* area_weights)

# end