import SimilarityWeights
using DimensionalData
using CairoMakie
using ColorSchemes
using NCDatasets

# load preprocessed data from ESMValTool #
#pathToPreprocDir = "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic"; 
pathToPreprocDir = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic"
climateVariables = ["pr"];
diagnostic = "";
weightsVarsIndep = Dict{String, Number}("pr" => 1); 
##########################################

modelData = SimilarityWeights.loadPreprocData(pathToPreprocDir, climateVariables, diagnostic, ["CMIP"]);
wI = SimilarityWeights.getIndependenceWeights(modelData, weightsVarsIndep);
wI

# plot wI
wI = dropdims(wI, dims=:variable);
f = SimilarityWeights.plotDistMatrices(wI, "pr", Array(dims(wI, :model1)), Array(dims(wI, :model2)))


pathObsData= "/albedo/work/projects/p_pool_clim_data/ERA5/Tier3/ERA5/v1/mon";
obsData = SimilarityWeights.loadPreprocData(pathObsData, ["tp"], diagnostic, ["era5"]);