import SimilarityWeights
using DimensionalData
using CairoMakie
using ColorSchemes
using NCDatasets

# load preprocessed data from ESMValTool #
pathToPreprocDir = "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic"; 
#pathToPreprocDir = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic"
climateVariables = ["pr"];
diagnostic = "";
weightsVars = Dict{String, Number}("pr" => 1); 
##########################################

modelData = SimilarityWeights.loadPreprocData(pathToPreprocDir, climateVariables, diagnostic, ["CMIP"]);
SimilarityWeights.averageEnsembleMembers!(modelData);
wI = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);
wI

# plot wI
wI = dropdims(wI, dims=:variable);
f = SimilarityWeights.plotDistMatrices(wI, "pr", Array(dims(wI, :model1)), Array(dims(wI, :model2)))

# on albedo
pathObsData= "/albedo/work/projects/p_pool_clim_data/ERA5/Tier3/ERA5/v1/mon";
obsData = SimilarityWeights.loadPreprocData(pathObsData, ["tp"], diagnostic, ["era5"]);

# local
pathObsData = "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240802_115940/preproc/climatologic_diagnostic";
obsData = SimilarityWeights.loadPreprocData(pathObsData, ["pr"], diagnostic, ["ERA5"]);
# just use one reference period
obsData["pr"] = obsData["pr"][:,:,1:1]; 

wP = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars)