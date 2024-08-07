import SimilarityWeights
using DimensionalData
using CairoMakie
using ColorSchemes
using NCDatasets

###########################################
pathsModelData = Dict("pr" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic");
pathsObsData = Dict("pr" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240802_115940/preproc/climatologic_diagnostic");
# on albedo
# pathObsData= "/albedo/work/projects/p_pool_clim_data/ERA5/Tier3/ERA5/v1/mon";
# obsData = SimilarityWeights.loadPreprocData(pathObsData, diagnostic, ["era5"]);
weightsVars = Dict{String, Number}("pr" => 1, "hur" => 1); 

#pathsModelData = Dict("hur" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_hur_20240801_094610/preproc/climatologic_diagnostic");
#pathsObsData = Dict("hur" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_hur_20240802_120722/preproc/climatologic_diagnostic");

## Get preprocessed data from ESMValTool ##
diagnostic = "";
modelData = SimilarityWeights.loadPreprocData(pathsModelData, diagnostic, ["CMIP"]);
obsData = SimilarityWeights.loadPreprocData(pathsObsData, diagnostic, ["ERA5"]);
# just use one reference period
obsData["pr"] = obsData["pr"][:,:,1:1]; 
###########################################

# Performance and independent weights
wP = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars);
wI = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);
wP_avg = SimilarityWeights.averageEnsembleVector(wP)
wI_avg = SimilarityWeights.averageEnsembleMatrix(wI)

# plot wI
wI = dropdims(wI_avg, dims=:variable);
f = SimilarityWeights.plotDistMatrices(wI, "pr", Array(dims(wI, :model1)), Array(dims(wI, :model2)))

weights = SimilarityWeights.getWeights(modelData, obsData)
