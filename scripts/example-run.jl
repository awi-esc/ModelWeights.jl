
using SimilarityWeights

path_config = "configs/example_historical_albedo.yml"
#path_config = "configs/example_historical_local.yml"
config = SimilarityWeights.validateConfig(path_config);

###############################################################################
#weights, avgs = SimilarityWeights.runWeights(config);
#SimilarityWeights.plotMeanData(config, avgs)

modelDataFull, modelDataRef, obsData = SimilarityWeights.getSharedModelData(config);

tos_data = modelDataFull["tos"]
data = SimilarityWeights.filterModels(tos_data, ["EC-Earth3"])
SimilarityWeights.plotEnsembleSpread(data, 7.5, 82.5)

###############################################################################



modelDataRef = SimilarityWeights.loadDataFromConfig(config, "name_ref_period", "models_project_name");
modelDataRef = SimilarityWeights.getCommonModelsAcrossVars(modelDataRef);
modelDataFull = SimilarityWeights.loadDataFromConfig(config, "name_full_period", "models_project_name");
modelDataFull = SimilarityWeights.getCommonModelsAcrossVars(modelDataFull);
obsData = SimilarityWeights.loadDataFromConfig(config, "name_ref_period", "obs_data_name");



# TODO: make sure that there is observational data for all variables for 
# the respective reference period, for now this is just assumed but in 
# some rare cases it may be wrong

# compute weights just for those models for which also historical reference period is available
shared_models = intersect(
    first(values(modelDataFull)).metadata["full_model_names"], 
    first(values(modelDataRef)).metadata["full_model_names"]
);
SimilarityWeights.keepModelSubset!(modelDataFull, shared_models);
SimilarityWeights.keepModelSubset!(modelDataRef, shared_models);
weights = SimilarityWeights.overallWeights(
    modelDataRef, 
    obsData, 
    config.weight_contributions["performance"],
    config.weight_contributions["independence"], 
    config.weights_variables["performance"],
    config.weights_variables["independence"]   
);   

means = SimilarityWeights.getWeightedAverages(modelDataFull, weights);
