import SimilarityWeights
using DimensionalData
using CairoMakie
using ColorSchemes
using NCDatasets

###########################################
pathsModelData = Dict("pr" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic");
pathsObsData = Dict("pr" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240802_115940/preproc/climatologic_diagnostic");

#pathsModelData = Dict("hur" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_hur_20240801_094610/preproc/climatologic_diagnostic");
#pathsObsData = Dict("hur" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_hur_20240802_120722/preproc/climatologic_diagnostic");

## Get preprocessed data from ESMValTool ##
#pathToPreprocDir = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic"
diagnostic = "";
modelData = SimilarityWeights.loadPreprocData(pathsModelData, diagnostic, ["CMIP"]);
# pathObsData = "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240802_115940/preproc/climatologic_diagnostic";
obsData = SimilarityWeights.loadPreprocData(pathsObsData, diagnostic, ["ERA5"]);
# just use one reference period
obsData["pr"] = obsData["pr"][:,:,1:1]; 
###########################################


weightsVars = Dict{String, Number}("pr" => 1, "hur" => 1); 
# Performance weights
wP = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars)

# Independence weights
SimilarityWeights.averageEnsembleMembers!(modelData);
wI = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);

# plot wI
wI = dropdims(wI, dims=:variable);
f = SimilarityWeights.plotDistMatrices(wI, "pr", Array(dims(wI, :model1)), Array(dims(wI, :model2)))

# on albedo
# pathObsData= "/albedo/work/projects/p_pool_clim_data/ERA5/Tier3/ERA5/v1/mon";
# obsData = SimilarityWeights.loadPreprocData(pathObsData, diagnostic, ["era5"]);


