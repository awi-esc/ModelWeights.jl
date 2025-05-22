@testset "Test buildCMIP5EnsembleMember" begin
    realizations = [1,2,3]
    initialization_methods = [2,3,4]
    physics = [1,1,1]
    member_ids = ModelWeights.buildCMIP5EnsembleMember(realizations, initialization_methods, physics)
    @test member_ids == ["r1i2p1", "r2i3p1", "r3i4p1"]
end

@testset "Test setLookupsFromMemberToModel" begin
end

@testset "Test reduceMetaDataSharedModels!" begin 
end

@testset "Test getModelsFromModelIDs" begin
end

@testset "Test getPhysicsFromMembers" begin
end