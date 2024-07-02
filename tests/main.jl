
import SimilarityWeights
using NCDatasets

# TODO: add real tests!

pathToPreprocDir = "recipe_climwip_test_basic_data/preproc/calculate_weights_climwip";

variables_long = ["Near-Surface Air Temperature", "Precipitation"];
climateVariables = ["tas", "pr"];
diagnostic = "CLIM";
included = ["CMIP"];
data = SimilarityWeights.loadPreprocData(pathToPreprocDir, climateVariables, diagnostic, included);



weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25); 


weights = SimilarityWeights.getIndependenceWeights(data, weightsVars);
Array(weights)

PATH_TO_WORK_DIR = joinpath(@__DIR__, "..", "recipe_climwip_test_basic_data", "work");

ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_overall_mean.nc"), "r");
nbDigits = 4;

dataTrue = round.(ds["overall_mean"][:,:], digits=nbDigits)
dataComputed = round.(Array(weights), digits=nbDigits)

@assert dataTrue == dataComputed
print("run main.jl")