@testset "Test combineAll" begin
    # unequal number of entries (first more)
    v1 = ["tos", "tas"]; v2 = ["CLIM"];
    expected1 = ["tas_CLIM", "tos_CLIM"]
    result1 = sort(ModelWeights.Data.combineAll(v1, v2))
    @test result1 == expected1

    # one/two input vectors empty -> should throw warning
    msg = "At least one input vector is empty -> empty vector returned!"
    @test_logs (:warn, msg) ModelWeights.Data.combineAll(Vector{String}(), Vector{String}())
    @test_logs (:warn, msg) ModelWeights.Data.combineAll(["a"], Vector{String}())

    # unequal number of arguments (second more)
    v1 = ["pr", "psl"]; v2 = ["CLIM", "std", "CLIM-mon"]
    result2 = sort(ModelWeights.Data.combineAll(v1, v2))
    expected2 = ["pr_CLIM", "pr_CLIM-mon", "pr_std", "psl_CLIM", "psl_CLIM-mon", "psl_std"]
    @test result2 == expected2
end


@testset "Test individuatePath" begin
end


@testset "Test renameDict" begin
    d = Dict{String, Vector{Int}}("b" => [7, 10], "c" => [17, 2])
    ModelWeights.Data.renameDict!(d, ["b"], ["bg"])
    @test haskey(d, "bg")
    @test !haskey(d, "b")

    d = Dict{Symbol, String}(:b => "oct", :c => "feb")
    ModelWeights.Data.renameDict!(d, [:c], [:cb])
    @test haskey(d, :cb)
    @test !haskey(d, :c)
end


@testset "Test kelvinToCelsius" begin
    mat = rand(2,3,4)
    mat[:,:, [1,3]] .+= 280
    mat[:,:, [2,4]] .+= 28 
    data = YAXArray(
        (Dim{:lon}([0, 100]), Dim{:lat}([-10, 0, 10]), Dim{:model}(["m1", "m2", "m3", "m4"])), 
        mat,
        Dict{String, Any}("units" => ["K", "degC", "K", "degC"])
    )
    df = ModelWeights.Data.kelvinToCelsius(data)
    @test all(df.properties["units"] .== "degC") # changes in new meta data
    @test data.properties["units"][[1,3]] == ["K", "K"] # no changes in original metadata
    @test all(df.data[:, :, [1, 3]] .< 280) # changes in new data
    @test all(data.data[:, :, [1, 3]] .> 280) # no changes in original data
end


@testset "Test dimNames" begin
    mat = rand(2,3,4)
    dimensions = (Dim{:lon}([0, 100]), Dim{:lat}([-10, 0, 10]), Dim{:model}(["m1", "m2", "m3", "m4"]))
    data = YAXArray(dimensions, mat)

    names = ModelWeights.Data.dimNames(data)
    @test all(x -> x in names, [:lon, :lat, :model])
    @test length(names) == 3
end

@testset "Test setDim sideeffects" begin
    mat = rand(2,3,4)
    dimensions = (Dim{:lon}([0, 100]), Dim{:lat}([-10, 0, 10]), Dim{:model}(["m1", "m2", "m3", "m4"]))
    
    data = YAXArray(dimensions, mat)
    df = ModelWeights.Data.setDim(data, :model, "mymodel", nothing)
    dimensions_data = ModelWeights.Data.dimNames(data)
    dimensions_df = ModelWeights.Data.dimNames(df)

    # test dimension names
    @test :model in dimensions_data
    @test !(:mymodel in dimensions_data)
    @test :mymodel in dimensions_df
    @test !(:model in dimensions_df)

    # test dimension values
    data = YAXArray(dimensions, mat)
    df = ModelWeights.Data.setDim(data, :model, "MODEL", ["M1", "M2", "M3", "M4"])
    @test lookup(df, :MODEL) == ["M1", "M2", "M3", "M4"]
    @test lookup(data, :model) == ["m1", "m2", "m3", "m4"]
end


@testset "Test countMap" begin
    data1 = ["a", "b", "b", "c", "b", "a"]
    counts = ModelWeights.Data.countMap(data1)
    @test counts["a"] == 2 && counts["b"] == 3 && counts["c"] == 1

    data2 = [7, 10, 17, 17, 91, 91, 910]
    counts = ModelWeights.Data.countMap(data2)
    @test counts[7] == 1 && counts[10] == 1 && counts[17] == 2 && counts[91] == 2 && counts[910] == 1
end