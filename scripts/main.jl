import SimilarityWeights
using DimensionalData
using CairoMakie
using ColorSchemes
using NCDatasets


RUN_ON_ALBEDO = true;

############# Get preprocessed data from ESMValTool ######################  
begin 
    if RUN_ON_ALBEDO
        # just use one reference period
        pathsObsData= Dict("pr" => "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_historical_pr_20240815_124334/preproc/climatology_historical1/pr");
        pathsModelData = Dict("pr" => "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_historical_pr_20240815_124334/preproc/climatology_historical1/pr");
    else
        pathsObsData = Dict("pr" => "/Users/brgrus001/output-from-albedo/recipe_historical_pr_20240815_124334/preproc/climatology_historical1/pr")
        pathsModelData = Dict("pr" => "/Users/brgrus001/output-from-albedo/recipe_historical_pr_20240815_124334/preproc/climatology_historical1/pr")
    end
end

obsData = SimilarityWeights.loadPreprocData(pathsObsData, ["ERA5"]);
weightsVars = Dict{String, Number}("pr" => 1, "tas" => 1, "hur" => 1); 
modelData = SimilarityWeights.loadPreprocData(pathsModelData, ["CMIP"]);
##########################################################################

# Performance and independent weights
wP = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars);
wI = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);
wP_avg = SimilarityWeights.averageEnsembleVector(wP, false)
wI_avg = SimilarityWeights.averageEnsembleMatrix(wI, false)

# plot wI
wI = dropdims(wI_avg, dims=:variable);
f = SimilarityWeights.plotDistMatrices(wI, "pr", Array(dims(wI, :model1)), Array(dims(wI, :model2)))

# plot wP
fig_wp = SimilarityWeights.plotPerformanceWeights(wP_avg)

weights = SimilarityWeights.getOverallWeights(modelData, obsData)
