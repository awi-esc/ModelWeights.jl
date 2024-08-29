include("data.jl")
using Statistics

@testset "Testset performance weights" begin    
    weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
    modelData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["CMIP"]);
    obsData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["ERA5"]);
    
    weights = SimilarityWeights.generalizedDistancesPerformance(modelData, obsData, weightsVars);
    weightsRounded = round.(Array(weights), digits=NB_DIGITS);

    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_overall_mean.nc"));
    expected =  Array(ds["overall_mean"]);
    expectedRounded = round.(expected, digits=NB_DIGITS);

    @test all((x)-> x==1, expectedRounded .== weightsRounded);
end


@testset "Testset independence weights" begin
    weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 
    modelData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["CMIP"]);

    weights = SimilarityWeights.generalizedDistancesIndependence(modelData, weightsVars);
    weights = round.(Array(weights), digits=NB_DIGITS);

    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_overall_mean.nc"), "r");
    expected = round.(ds["overall_mean"][:,:], digits=NB_DIGITS);

    @test all((x)-> x==1, expected .== weights);
end


# TODO: go through these tests and adapt!


# @testset "Test averaging ensemble members Matrix (independence weights)" begin
#     modelData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["CMIP"]);
#     weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 

#     weights = SimilarityWeights.generalizedDistancesIndependence(modelData, weightsVars);
#     wInd_vars = SimilarityWeights.averageEnsembleMatrix(weights)
#     wInd = SimilarityWeights.summarizeWeightsAcrossVars(wInd_vars)

#     expected_all = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_overall_mean.nc"))["overall_mean"];    
#     #last four entries are from one ensemble 
#     expected = copy(expected_all[1:4, 1:4])
#     expected[4,4] = 0
#     for i in 1:3
#         val = mean(expected_all[i, 4:end]);
#         expected[i, 4] = val;
#         expected[4, i] = val;
#     end
#     @test all((x)-> x==1, round.(expected, digits=NB_DIGITS) .== round.(Array(wInd), digits=NB_DIGITS));

# end

# @testset "Test averaging ensemble members Vector (performance weights)" begin
#     modelData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["CMIP"]);
#     weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 

#     obsData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["ERA5"]);
#     weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
#     weights = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars);
#     wPerform_vars = SimilarityWeights.averageEnsembleVector(weights);
#     wPerform = SimilarityWeights.summarizeWeightsAcrossVars(wPerform_vars);
    
#     expected_all = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_overall_mean.nc"))["overall_mean"];    
#     #last four entries are from one ensemble 
#     expected = copy(expected_all[1:3])
#     push!(expected, mean(expected_all[4:end]))
#     @test all((x)-> x==1, round.(expected, digits=NB_DIGITS) .== round.(Array(wPerform), digits=NB_DIGITS));
# end


# @testset "Testset combined weights" begin
#     modelData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["CMIP"]);
#     weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 
#     independenceWeights = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);

#     obsData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["ERA5"]);
#     weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
#     performanceWeights = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars);
#     weights = SimilarityWeights.combineWeights(performanceWeights, independenceWeights, 0.5, 0.5)

#     ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "weights.nc"));
#     expected_all = Array(ds["weight"]);
#     #last four entries are from one ensemble 
#     expected = copy(expected_all[1:3])
#     push!(expected, sum(expected_all[4:end]))

#     @test all((x)-> x==1, round.(expected, digits=NB_DIGITS) .== round.(weights, digits=NB_DIGITS));
# end


@testset "Testset function getOverallWeights" begin
    modelData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["CMIP"]);
    obsData = SimilarityWeights.loadPreprocData(VAR_TO_PREPROC_DATA, ["ERA5"]);
    
    weightsVarsPerform = Dict{String, Number}("tas" => 1, 
                                              "pr" => 2, 
                                              "psl" => 1); 
    weightsVarsIndep = Dict{String, Number}("tas" => 0.5, 
                                            "pr" => 0.25,
                                            "psl" => 0); 
    weights = SimilarityWeights.overallWeights(
        modelData, obsData, 0.5, 0.5, weightsVarsPerform, weightsVarsIndep
    );

    ds = NCDataset(joinpath(PATH_TO_WORK_DIR,
                   "calculate_weights_climwip", 
                   "climwip", 
                   "weights.nc"));
    expected_all = Array(ds["weight"]);

    #last four entries are from one ensemble 
    expected = copy(expected_all[1:3]);
    push!(expected, sum(expected_all[4:end]))
 
    @test all((x)-> x==1, round.(expected, digits=NB_DIGITS) .== round.(weights, digits=NB_DIGITS));
end


