@testset "Test buildCMIP5EnsembleMember" begin
    realizations = [1,2,3]
    initialization_methods = [2,3,4]
    physics = [1,1,1]
    vals = collect(zip(realizations, initialization_methods, physics))
    member_ids = map(x -> ModelWeights.Data.buildCMIP5EnsembleMember.(x...), vals)
    @test member_ids == ["r1i2p1", "r2i3p1", "r3i4p1"]
end

@testset "Test setLookupsFromMemberToModel" begin
end

@testset "Test reduceMetaDataSharedModels!" begin 
end

@testset "Test getModelsFromModelIDs" begin
end

@testset "Test physicsFromMember" begin
end