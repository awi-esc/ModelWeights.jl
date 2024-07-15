include("data.jl")

# CLIMATE_VARS = ["tas"];
# DIAGNOSTIC = "CLIM";
# MODEL_DATA =  SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, CLIMATE_VARS, DIAGNOSTIC, ["CMIP"]);
# OBS_DATA = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, CLIMATE_VARS, DIAGNOSTIC, ["ERA5"]);

# VARIABLE_CONTRIBUTIONS = Dict(
#     "independence" => Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0), 
#     "performance" => Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1)
# );

@testset "Testset performance weights" begin    
    # weights = SimilarityWeights.getPerformanceWeights(MODEL_DATA, OBS_DATA, VARIABLE_CONTRIBUTIONS["performance"]);
    weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
    modelData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["CMIP"]);
    obsData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["ERA5"]);
    
    weights = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars);
    weightsRounded = round.(Array(weights), digits=NB_DIGITS);

    #weightsEnsemble = SimilarityWeights.averageEnsembleVector(weights);
    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_overall_mean.nc"));
    expected =  Array(ds["overall_mean"]);
    expectedRounded = round.(expected, digits=NB_DIGITS);

    @test all((x)-> x==1, expectedRounded .== weightsRounded);
end


@testset "Testset independence weights" begin
    # weights = SimilarityWeights.getIndependenceWeights(MODEL_DATA, VARIABLE_CONTRIBUTIONS["independence"]);
    weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 
    modelData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["CMIP"]);
    weights = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);
    weights = round.(Array(weights), digits=NB_DIGITS);

    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_overall_mean.nc"), "r");
    expected = round.(ds["overall_mean"][:,:], digits=NB_DIGITS);

    @test all((x)-> x==1, expected .== weights);
end


# @testset "Testset combined weights" begin
#     # weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 
#     wI = SimilarityWeights.getIndependenceWeights(MODEL_DATA, VARIABLE_CONTRIBUTIONS["independence"]);

#     # weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
#     wP = SimilarityWeights.getPerformanceWeights(MODEL_DATA, OBS_DATA, VARIABLE_CONTRIBUTIONS["performance"]);
#     weights = SimilarityWeights.combineWeights(wP, wI, 0.5, 0.5)

#     ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "weights.nc"));
#     expected = ds["weight"];

#     @test round.(weights, digits=NB_DIGITS) == round.(expected, digits=NB_DIGITS)
# end