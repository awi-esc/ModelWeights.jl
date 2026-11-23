import DimensionalData as DD
import ModelWeights.Data as mwd

using DimensionalData: Dim, dims, lookup
using Test
using YAXArrays

const CLIMWIP_DATA = mwd.loadDataFromESMValToolRecipes(
    joinpath(@__DIR__, "data/climwip-simplified"),
    joinpath(@__DIR__, "../configs/climwip_config");
    recipe_per_var = false,
    filename_format = :esmvaltool,
    constraint = Dict(
        :aliases => ["calculate_weights_climwip"],
        :mips => ["CMIP5", "CMIP6"],
    ),
)

@testset "Test combineAll" begin
    # unequal number of entries (first more)
    v1 = ["tos", "tas"]
    v2 = ["CLIM"]
    expected1 = ["tas_CLIM", "tos_CLIM"]
    result1 = sort(mwd.combineAll(v1, v2))
    @test result1 == expected1

    # one/two input vectors empty -> should throw warning
    msg = "At least one input vector is empty -> empty vector returned!"
    @test_logs (:warn, msg) mwd.combineAll(Vector{String}(), Vector{String}())
    @test_logs (:warn, msg) mwd.combineAll(["a"], Vector{String}())

    # unequal number of arguments (second more)
    v1 = ["pr", "psl"]
    v2 = ["CLIM", "std", "CLIM-mon"]
    result2 = sort(mwd.combineAll(v1, v2))
    expected2 = ["pr_CLIM", "pr_CLIM-mon", "pr_std", "psl_CLIM", "psl_CLIM-mon", "psl_std"]
    @test result2 == expected2
end


# @testset "Test individuatePath" begin
# end


@testset "Test renameDict" begin
    d = Dict{String,Vector{Int}}("b" => [7, 10], "c" => [17, 2])
    mwd.renameDict!(d, ["b"], ["bg"])
    @test haskey(d, "bg")
    @test !haskey(d, "b")

    d = Dict{Symbol,String}(:b => "oct", :c => "feb")
    mwd.renameDict!(d, [:c], [:cb])
    @test haskey(d, :cb)
    @test !haskey(d, :c)
end


@testset "Test kelvinToCelsius" begin
    mat = rand(2, 3, 4)
    mat[:, :, [1, 3]] .+= 280
    mat[:, :, [2, 4]] .+= 28
    data = YAXArray(
        (
            Dim{:lon}([0, 100]),
            Dim{:lat}([-10, 0, 10]),
            Dim{:model}(["m1", "m2", "m3", "m4"]),
        ),
        mat,
        Dict{String,Any}("units" => ["K", "degC", "K", "degC"]),
    )
    df = mwd.kelvinToCelsius(data)
    @test all(df.properties["units"] .== "degC") # changes in new meta data
    @test data.properties["units"][[1, 3]] == ["K", "K"] # no changes in original metadata
    @test all(df.data[:, :, [1, 3]] .< 280) # changes in new data
    @test all(data.data[:, :, [1, 3]] .> 280) # no changes in original data
end


@testset "Test dimNames" begin
    mat = rand(2, 3, 4)
    dimensions = (
        Dim{:lon}([0, 100]),
        Dim{:lat}([-10, 0, 10]),
        Dim{:model}(["m1", "m2", "m3", "m4"]),
    )
    data = YAXArray(dimensions, mat)

    names = mwd.dimNames(data)
    @test all(x -> x in names, [:lon, :lat, :model])
    @test length(names) == 3
end

@testset "Test setDim sideeffects" begin
    mat = rand(2, 3, 4)
    dimensions = (
        Dim{:lon}([0, 100]),
        Dim{:lat}([-10, 0, 10]),
        Dim{:model}(["m1", "m2", "m3", "m4"]),
    )

    data = YAXArray(dimensions, mat)
    df = mwd.setDim(data, :model, :mymodel)
    dimensions_data = mwd.dimNames(data)
    dimensions_df = mwd.dimNames(df)

    # test dimension names
    @test :model in dimensions_data
    @test !(:mymodel in dimensions_data)
    @test :mymodel in dimensions_df
    @test !(:model in dimensions_df)

    # test dimension values
    data = YAXArray(dimensions, mat)
    df = mwd.setDim(data, :model, :MODEL, ["M1", "M2", "M3", "M4"])
    @test lookup(df, :MODEL) == ["M1", "M2", "M3", "M4"]
    @test lookup(data, :model) == ["m1", "m2", "m3", "m4"]
end


@testset "Test countMap" begin
    data1 = ["a", "b", "b", "c", "b", "a"]
    counts = mwd.countMap(data1)
    @test counts["a"] == 2 && counts["b"] == 3 && counts["c"] == 1

    data2 = [7, 10, 17, 17, 91, 91, 910]
    counts = mwd.countMap(data2)
    @test counts[7] == 1 &&
          counts[10] == 1 &&
          counts[17] == 2 &&
          counts[91] == 2 &&
          counts[910] == 1
end


@testset "Test mergeYAX new dimension" begin
    longitudes = [12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5]
    latitudes = [-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5]
    dimensions = (
        Dim{:lat}(latitudes),
        Dim{:lon}(longitudes),
        Dim{:member}(["ESM1#r1i1p1f1", "ESM2#r1i1p1f1", "ESM2#r1i1p1f2"]),
    )
    arr1 = YAXArray(dimensions, zeros(7, 9, 3))
    arr2 = YAXArray(dimensions, ones(7, 9, 3))
    merged = mwd.mergeYAX([arr1, arr2], :var, ["tos", "tas"])
    @test size(merged) == (7, 9, 3, 2)
    dim_new = dims(merged, 4)
    @test DD.name(dim_new) == :var
    @test collect(dim_new) == ["tas", "tos"] # is always added as ForwardOrdered
end


@testset "Test mergeYAX extend dimension sorted" begin
    longitudes = [12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5]
    latitudes = [-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5]
    dims1 =
        (Dim{:lat}(latitudes), Dim{:lon}(longitudes), Dim{:model}(["ESM1", "ESM2", "ESM5"]))
    dims2 = (Dim{:lat}(latitudes), Dim{:lon}(longitudes), Dim{:model}(["ESM3", "ESM4"]))
    arr1 = YAXArray(dims1, zeros(7, 9, 3))
    arr2 = YAXArray(dims2, ones(7, 9, 2))
    df = mwd.mergeYAX(arr1, arr2, :model)

    @test sum(df[model=At("ESM3")]) == 7*9
    @test sum(df[model=3]) == 7*9

    @test sum(df[model=At("ESM5")]) == 0
    @test sum(df[model=5]) == 0

    @test lookup(df, :model) == ["ESM1", "ESM2", "ESM3", "ESM4", "ESM5"]

    # TODO: test if ForwardOrdered
end


@testset "Test mergeYAX extend dimension unsorted" begin
    longitudes = [12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5]
    latitudes = [-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5]
    dims1 =
        (Dim{:lat}(latitudes), Dim{:lon}(longitudes), Dim{:model}(["ESM1", "ESM2", "ESM5"]))
    dims2 = (Dim{:lat}(latitudes), Dim{:lon}(longitudes), Dim{:model}(["ESM3", "ESM4"]))
    arr1 = YAXArray(dims1, zeros(7, 9, 3))
    arr2 = YAXArray(dims2, ones(7, 9, 2))
    df = mwd.mergeYAX(arr1, arr2, :model; sorted = false)

    @test sum(df[model=At("ESM3")]) == 7*9
    @test sum(df[model=4]) == 7*9

    @test sum(df[model=At("ESM5")]) == 0
    @test sum(df[model=3]) == 0

    @test lookup(df, :model) == ["ESM1", "ESM2", "ESM5", "ESM3", "ESM4"]

    # TODO: test if UnOrdered
end

@testset "Test subsetModelData(::YAXArray)" begin
    df_models =
        mwd.summarizeMembers(CLIMWIP_DATA.models["pr_CLIM_calculate_weights_climwip"])
    df_members = CLIMWIP_DATA.models["pr_CLIM_calculate_weights_climwip"]
    n_models = length(dims(df_models, :model))

    df_sub = mwd.subsetModelData(df_models, ["CCSM4", "ACCESS1-0#r1i1p1"])
    @test length(dims(df_sub, :model)) == 1

    # data on level of models, filter on member level
    df_sub = mwd.subsetModelData(df_models, ["ACCESS1-0#r1i1p1"])
    @test ndims(df_sub) == 2

    df_sub = mwd.subsetModelData(df_members, ["CCSM4"])
    @test length(dims(df_sub, :member)) == 4

    df_sub = mwd.subsetModelData(df_members, ["CCSM4#r3i1p1"])
    @test length(dims(df_sub, :member)) == 1

    # subset specific run of a model and that model generally -> all runs should be kept
    df_sub = mwd.subsetModelData(df_members, ["CCSM4#r3i1p1", "CCSM4"])
    @test length(dims(df_sub, :member)) == 4

    df_sub = mwd.subsetModelData(df_models, ["AWI-ESM"])
    @test ndims(df_sub) == 2
end


@testset "Test subsetModelData(::DataMap)" begin
    dm = CLIMWIP_DATA.models
    df = dm["pr_CLIM_calculate_weights_climwip"]
    models = collect(lookup(df, :member))[1:(end-2)]

    dm["pr_CLIM_calculate_weights_climwip"] = mwd.subsetModelData(df, models)

    # nothing should change
    ids = collect(keys(dm))
    ids_subset = ["tas_CLIM_calculate_weights_climwip", "pr_CLIM_calculate_weights_climwip"]
    dm_sub = mwd.subsetModelData(dm; level = :model, ids = ids_subset)
    @test length(collect(keys(dm_sub))) == 2
    for id in ids_subset
        @test lookup(dm_sub[id], :member) == lookup(dm[id], :member)
    end
    # use a non-existent id
    @test_throws ArgumentError mwd.subsetModelData(dm; level = :model, ids = ["pr"])

    # check that nothing changes when using a single id
    dm_sub = mwd.subsetModelData(
        dm;
        level = :member,
        ids = ["pr_CLIM_calculate_weights_climwip"],
    )
    @test lookup(dm_sub["pr_CLIM_calculate_weights_climwip"], :member) ==
          lookup(dm["pr_CLIM_calculate_weights_climwip"], :member)

    # when subsetting on level :member, from the CCSM4 models, only two should remain everywhere
    dm_sub = mwd.subsetModelData(dm; level = :member, ids = ids_subset)
    @test length(collect(keys(dm_sub))) == 2
    for id in ids_subset
        @test size(dm_sub[id], :member) == 5
    end
    # and without using ids
    dm_sub = mwd.subsetModelData(dm; level = :member)
    @test length(collect(keys(dm_sub))) == 3
    for id in ids
        @test size(dm_sub[id], :member) == 5
    end
end
