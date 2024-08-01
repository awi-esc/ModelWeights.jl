import SimilarityWeights
using DimensionalData
using CairoMakie
using ColorSchemes

# load preprocessed data from ESMValTool #
pathToPreprocDir = "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic"; 
climateVariables = ["pr"];
diagnostic = "";
weightsVarsIndep = Dict{String, Number}("pr" => 1); 
##########################################


modelData = SimilarityWeights.loadPreprocData(pathToPreprocDir, climateVariables, diagnostic, ["CMIP"]);
wI = SimilarityWeights.getIndependenceWeights(modelData, weightsVarsIndep);
wI