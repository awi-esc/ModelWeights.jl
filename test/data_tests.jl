import ModelWeights.Data as mwd

using DimensionalData: Dim, dims
using Test
using YAXArrays


longitudes = [12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5]
latitudes = [-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5]
arr1 = YAXArray(
    (
        Dim{:lat}(latitudes),
        Dim{:lon}(longitudes),
        Dim{:member}(["ESM1#r1i1p1f1", "ESM2#r1i1p1f1", "ESM2#r1i1p1f2"]),
    ),
    rand(7, 9, 3),
)
arr2 = YAXArray(
    (
        Dim{:lat}(latitudes),
        Dim{:lon}(longitudes),
        Dim{:member}(["ESM1#r1i1p1f1", "ESM2#r1i1p1f1"]),
    ),
    rand(7, 9, 2),
)

weights = YAXArray(
    (Dim{:member}(["ESM1#r1i1p1f1", "ESM2#r1i1p1f1", "ESM2#r1i1p1f2"]),),
    [0.2, 0.3, 0.5],
)

a = [1.0 2.0 missing 4];
b = [5.0 missing missing 8];
c = [missing -1 0 1];

d1 = reshape(vcat([1 2], [3 4], [5 6]), 2, 3, 1);
d2 = reshape(1:6, 2, 3, 1);
yax1 = YAXArray(
    (Dim{:row}(["r1", "r2"]), Dim{:column}(["c1", "c2", "c3"]), Dim{:stack}(["s1"])),
    d1,
)
yax2 = YAXArray(
    (Dim{:row}(["r1", "r2"]), Dim{:column}(["c1", "c2", "c3"]), Dim{:stack}(["s1"])),
    d2,
)


function setupDataPaths()
    paths_lgm_tas = ["data/lgm-cmip5-tas-climatologies", "data/lgm-cmip6-tas-climatologies"]
    paths_lgm_tos = ["data/lgm-cmip5-tos-climatologies", "data/lgm-cmip6-tos-climatologies"]
    paths_lgm_tas = joinpath.(pwd(), paths_lgm_tas)
    paths_lgm_tos = joinpath.(pwd(), paths_lgm_tos)
    return Dict("lgm-tas" => paths_lgm_tas, "lgm-tos" => paths_lgm_tos)
end


@testset "Test joinDataMaps" begin
    dm1 = mwd.DataMap(Dict("id1" => yax1))
    dm2 = mwd.DataMap(Dict("id2" => yax2))
    @test mwd.joinDataMaps(dm1, dm2) == mwd.DataMap(Dict("id1" => yax1, "id2" => yax2))
    # warning is thrown and value of shared key is taken from the second argument
    @test (@test_logs (:warn,) mwd.joinDataMaps(
        dm1,
        mwd.DataMap(Dict("id1" => yax2));
        warn_msg = "my warning",
    )) == mwd.DataMap(Dict("id1" => yax2))
end

@testset "Test metaDataFromESMValToolRecipes" begin end

@testset "Test metaDataFromYAML" begin end

@testset "Test resolvePathsFromMetaData" begin end


@testset "Test loadDataFromYAML" begin
    # path = "/albedo/home/brgrus001/ModelWeightsPaper/work/configs/climwip/climwip-model-data.yml"
    # data = mwd.loadDataFromYAML(path, preview = true);
end


@testset "Test getAtModel" begin
    a = [1 2 3]
    b = [4 5 6]
    da = YAXArray((Dim{:model}(["m1", "m2"]), Dim{:var}(["tas", "tos", "pr"])), vcat(a, b))
    @test mwd.getAtModel(da, "m1") == [1, 2, 3]

    mwd.putAtModel!(da, "m2", [17, 2, 1987])
    @test mwd.getAtModel(da, "m2") == [17, 2, 1987]
end

# TODO: The following testset retrieves local data, update this to use on-the-fly generated pseudo data.
# @testset "Test defineDataMap(path::String, id::String) filename_format :cmip" begin
#     # Original (from ESGF) CMIP data, may include the timerange or not, for both cmip5 and cmip6.

#     # directory contains pr and tas data for CCSM4 model. The distinguishing variable is not model, but variable...
#     # TODO: this is not supported! (does run, but loads metadata for both variables, but just the first variable's data)
#     path_cmip5 = "/albedo/work/projects/p_pool_clim_data/cmip5/output1/NCAR/CCSM4/midHolocene/mon/atmos/Amon/r1i1p1/v20120604";

#     for fn in [:cmip, :cmip5]
#         df = mwd.defineDataMap(
#             path_cmip5,
#             "midHolocene";
#             filename_format = :cmip,
#             constraint = Dict(:variables => ["tas"])
#         )
#         @test Array(df["midHolocene"].member) == ["CCSM4#r1i1p1"]
#         @test df["midHolocene"].properties["name"] == ["tas"]
#     end

#     # do not constrain variables:
#     df = mwd.defineDataMap(
#         path_cmip5,
#         "midHolocene";
#         filename_format = :cmip
#     )
#     # since the directory contains multiple variables for the same model, it should return all data,
#     # and model CCSM4 should be twice in the model dimension. TODO: this is not supported yet
#     @test_broken length(dims(df["midHolocene"], :member)) == 2
#     @warn "The directory should contain data for one variable and multiple models, or just one model and variable for multiple time periods, but not multiple variables for one model."

#     # CMIP6
#     # directory only contains orog_fx_AWI-ESM-1-REcoM_historical_r1i1p1f1_gn.nc
#     path_cmip6 = "/albedo/work/projects/p_pool_clim_data/CMIP6/CMIP/AWI/AWI-ESM-1-REcoM/historical/r1i1p1f1/fx/orog/gn/v20230314/"
#     for fn in [:cmip, :cmip6]
#         df = mwd.defineDataMap(
#             path_cmip6,
#             "orog-historical";
#             filename_format = fn
#         )
#         @test Array(df["orog-historical"].member) == ["AWI-ESM-1-REcoM#r1i1p1f1"]
#     end
#     # directory contains multiple files with pr for AWI-ESM-1-REcoM, for various time periods, from 1850-2014
#     path_cmip6 = "/albedo/work/projects/p_pool_clim_data/CMIP6/CMIP/AWI/AWI-ESM-1-REcoM/historical/r1i1p1f1/Amon/pr/gn/v20230314"
#     for fn in [:cmip, :cmip6_tr]
#         df = mwd.defineDataMap(
#             path_cmip6,
#             "pr-historical-ts";
#             filename_format = fn
#         )
#         @test year(df["pr-historical-ts"].time[1]) == 1850
#         @test year(df["pr-historical-ts"].time[end]) == 2014
#     end
# end

# Needs to be updated:
# @testset "Test defineDataMap" begin
#     filename_format = :esmvaltool
#     dtype = "cmip"
#     data_paths = setupDataPaths()

#     paths_lgm = [data_paths["lgm-tas"], data_paths["lgm-tos"]]
#     ids = ["lgm-tas-data", "tos_lgm"]
#     data = mwd.defineDataMap(paths_lgm, ids; dtype, filename_format)
#     @test size(data["lgm-tas-data"]) == (72, 36, 17)
#     @test size(data["tos_lgm"]) == (72, 36, 16)

#     constraint = Dict{String, Union{Vector{String}, Symbol}}("level_shared" => :member)
#     data = mwd.defineDataMap(paths_lgm, ids; dtype, filename_format, constraint)
#     @test data["lgm-tas-data"].member == data["tos_lgm"].member
#     @test length(data["lgm-tas-data"].member) == 15

#     constraint["models"] = ["GISS-E2-R"]
#     data = mwd.defineDataMap(paths_lgm, ids; dtype, filename_format, constraint)
#     @test length(data["lgm-tas-data"].member) == 2

#     constraint["models"] = ["GISS-E2-R#r1i1p150"]
#     data = mwd.defineDataMap(paths_lgm, ids; dtype, filename_format, constraint)
#     @test length(data["lgm-tas-data"].member) == 1
# end
