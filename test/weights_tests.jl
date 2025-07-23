import ModelWeights.Data as mwd
import ModelWeights.Weights as mww

@testset "Test distancesData" begin
end

@testset "Test distancesModels" begin
end

@testset "Test generalizedDistances" begin
end

@testset "Test uncertaintyRanges" begin
end

@testset "Test climwipWeights" begin
end

@testset "Test areaWeightedRMSE" begin
end

@testset "Test normalize" begin
end

@testset "Test summarizeMembers vectors" begin
    members = ["ACCESS-CM2#r10i1p1f1_gn", "CNRM-CM6-1#r13i1p1f2_gr", "CNRM-CM6-1#r15i1p1f2_gr", "MIROC6#r23i1p1f1_gn"]
    vals = [1, 2, 3, 4, 5.5, 6.5, 7.5, 8.5]
    meta = Dict(
        "id" => ["ACCESS-CM2", "CNRM-CM6-1#r13i1p1f2_gr", "CNRM-CM6-1#r15i1p1f2_gr", "MIROC6"],
        "variable" => ["tas", "tas", "tas", "tas"],
        "info" => Dict("any-key" => "any-value")
    )
    data = YAXArray((
        Dim{:member}(members), Dim{:time}(["Jan", "Feb"])), reshape(vals, 4, 2), meta
    )
    averaged = mwd.summarizeMembers(data)
    @test Array(averaged[model = At("CNRM-CM6-1")]) == [2.5, 7.0]
    @test length(averaged.model) == 3
    # test if metadata also has the correct size
    @test length(averaged.properties["id"]) == 3
    @test length(averaged.properties["id"][2]) == 2
    @test length(averaged.properties["variable"]) == 3
    
    info = get(averaged.properties, "info", nothing)
    @test !isnothing(info) && !isempty(info)
end

@testset "Test summarizeMembers matrices" begin
end

@testset "Test weightedAvg" begin
end

@testset "Test equalWeights" begin
    members = ["ACCESS-CM2#r10i1p1f1_gn", "CNRM-CM6-1#r13i1p1f2_gr", "CNRM-CM6-1#r15i1p1f2_gr", "MIROC6#r23i1p1f1_gn"]
    vals = [1, 2, 3, 4, 5.5, 6.5, 7.5, 8.5]
    data = YAXArray((Dim{:member}(members), Dim{:time}(["Jan", "Feb"])),reshape(vals, 4, 2))
    w1 = mww.equalWeights(data;use_members=true)
    @test Array(w1) == [1/3, 1/6, 1/6, 1/3]
    w2 = mww.equalWeights(data;use_members=false)
    @test Array(w2) == [1/4, 1/4, 1/4, 1/4]

    models = ["ACCESS-CM2", "CNRM-CM6-1", "MIROC6"]
    vals = [1, 2, 3, 4, 5, 6]
    data = YAXArray(
        (Dim{:model}(models), Dim{:time}(["Jan", "Feb"])),
        reshape(vals, 3, 2)
    )
    w = mww.equalWeights(data)
    @test Array(w) == [1/3, 1/3, 1/3]
end

@testset "Test distributeWeightsAcrossMembers" begin
    models = ["ACCESS-CM2", "CNRM-CM6-1", "MIROC6"]
    members = ["ACCESS-CM2#r10i1p1f1_gn", "CNRM-CM6-1#r13i1p1f2_gr", 
        "CNRM-CM6-1#r15i1p1f2_gr", "MIROC6#r23i1p1f1_gn"]
    w = YAXArray((Dim{:model}(models),),[0.1, 0.5, 0.4])
    w_members = mww.distributeWeightsAcrossMembers(w, members)
    @test w_members == [0.1, 0.25, 0.25, 0.4]
end

@testset "Test saveWeightsAsNCFile" begin  
end

@testset "Test saveWeightsAsJuliaObj" begin 
end

@testset "Test loadWeightsFromJLD2" begin
end

@testset "Test applyWeights" begin 
end