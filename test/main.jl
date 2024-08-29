import SimilarityWeights
using NCDatasets
using DimensionalData
using Statistics


# TODO: go through this script not up to date anymore!

pathToPreprocDir = "recipe_climwip_test_basic_data/preproc/calculate_weights_climwip";
PATH_TO_WORK_DIR = joinpath(@__DIR__, "..", "recipe_climwip_test_basic_data", "work");

climateVariables = ["tas", "pr", "psl"];
diagnostic = "CLIM";

modelData = SimilarityWeights.loadPreprocData(pathToPreprocDir, climateVariables, diagnostic, ["CMIP"]);
obsData = SimilarityWeights.loadPreprocData(pathToPreprocDir, climateVariables, diagnostic, ["ERA5"]);

# 1. Independence weights
weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 
independenceWeights = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);

# test
nbDigits = 6;

ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_overall_mean.nc"), "r");
independenceTrue = round.(ds["overall_mean"][:,:], digits=nbDigits);
independenceComputed = round.(Array(independenceWeights), digits=nbDigits);

independenceTrue .== independenceComputed

@assert !any((x) -> x!=1, independenceTrue .== independenceComputed)
# ########################################
# # 2. Performance weights
climVar = "pr";

performanceDists = SimilarityWeights.getModelDataDist(modelData[climVar], obsData[climVar]);
# compare them to true data (from work dir)
ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_" * climVar * "_CLIM.nc"));
performanceDistTrue = ds["d" * climVar * "_CLIM"];

round.(Array(performanceDistTrue), digits=nbDigits) .== round.(Array(performanceDists), digits = nbDigits)

@assert round.(Array(performanceDistTrue), digits=nbDigits) == round.(Array(performanceDists), digits=nbDigits);


weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
performanceWeights = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars);
performanceComputed = round.(Array(performanceWeights), digits=nbDigits);

ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_overall_mean.nc"));
performanceTrue = round.(Array(ds["overall_mean"]), digits=nbDigits);

performanceTrueByEnsemble = performanceTrue[1:4];
performanceTrueByEnsemble[4] =  Statistics.mean(performanceTrue[4:end]);
performanceTrueByEnsemble = round.(performanceTrueByEnsemble, digits = nbDigits);

@assert !any((x) -> x!=1, performanceTrue .== performanceComputed)


# #### Combine weights
weights = SimilarityWeights.combineWeights(performanceWeights, independenceWeights, 0.5, 0.5)

println("finished")