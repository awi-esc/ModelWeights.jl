import SimilarityWeights
using DimensionalData
using CairoMakie
using ColorSchemes
using NCDatasets


RUN_ON_ALBEDO = true;

############# Get preprocessed data from ESMValTool ######################
if RUN_ON_ALBEDO
    # just use one reference period
    pathsObsData= Dict("pr" => "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_historical_pr_20240815_124334/preproc/climatology_historical1/pr");
    pathsModelData = Dict("pr" => "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_historical_pr_20240815_124334/preproc/climatology_full/pr");
else
    # TODO: recheck these!!
    pathsObsData = Dict("tp" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240802_115940/preproc/climatologic_diagnostic",
                        "hur" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_hur_20240802_120722/preproc/climatologic_diagnostic");
    pathsModelData = Dict("pr" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic", 
                          "hur" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_hur_20240801_094610/preproc/climatologic_diagnostic");
end

obsData = SimilarityWeights.loadPreprocData(pathsObsData, ["ERA5"]);
weightsVars = Dict{String, Number}("pr" => 1, "tas" => 1, "hur" => 1); 
modelData = SimilarityWeights.loadPreprocData(pathsModelData, ["CMIP"]);
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
