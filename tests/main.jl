
import SimilarityWeights
using NCDatasets
using DimensionalData
using Statistics
# TODO: add real tests!

pathToPreprocDir = "recipe_climwip_test_basic_data/preproc/calculate_weights_climwip";
PATH_TO_WORK_DIR = joinpath(@__DIR__, "..", "recipe_climwip_test_basic_data", "work");

climateVariables = ["tas", "pr", "psl"];
diagnostic = "CLIM";
included = ["CMIP"];
weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 

modelData = SimilarityWeights.loadPreprocData(pathToPreprocDir, climateVariables, diagnostic, included);

# 1. Independence weights
weights = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);
Array(weights)

# test
nbDigits = 5;

ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_overall_mean.nc"), "r");
independenceTrue = round.(ds["overall_mean"][:,:], digits=nbDigits);
independenceComputed = round.(Array(weights), digits=nbDigits);

independenceTrue .== independenceComputed

@assert !any((x) -> x!=1, independenceTrue .== independenceComputed)
########################################
# 2. Performance weights
obsData = SimilarityWeights.loadPreprocData(pathToPreprocDir, climateVariables, diagnostic, ["ERA5"]);
# Make sure to use a copy of the observational data, otherwise, it will be modified within the function by applying the mask!!
climVar = climateVariables[2];
weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 


models = modelData[climVar];
observations = deepcopy(obsData[climVar]);
performanceDists = SimilarityWeights.getModelDataDist(models, observations);
# compare them to data from work dir
ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_" * climVar * "_CLIM.nc"));
performanceDistTrue = ds["d" * climVar * "_CLIM"];

round.(Array(performanceDistTrue), digits=nbDigits) .== round.(Array(performanceDists), digits = nbDigits)

@assert round.(Array(performanceDistTrue), digits=nbDigits) == round.(Array(performanceDists), digits=nbDigits);
#@assert Array(performanceDistTrue) == Array(performanceDists);


performanceWeights = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars);
performanceComputed = round.(Array(performanceWeights), digits=nbDigits);

ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_overall_mean.nc"));
performanceTrue = Array(ds["overall_mean"]);
performanceTrueByEnsemble = performanceTrue[1:4];
performanceTrueByEnsemble[4] =  Statistics.mean(performanceTrue[4:end]);
performanceTrue = round.(performanceTrueByEnsemble, digits = nbDigits);

@assert !any((x) -> x!=1, performanceTrue .== performanceComputed)

println("finished")