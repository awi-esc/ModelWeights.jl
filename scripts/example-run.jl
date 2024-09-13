
using SimilarityWeights

path_config = "configs/example_historical_albedo.yml"
path_config = "configs/example_historical_local.yml"
config = SimilarityWeights.validateConfig(path_config);

###############################################################################
#weights, avgs = SimilarityWeights.runWeights(config);
#SimilarityWeights.plotMeanData(config, avgs)

modelDataFull, modelDataRef, obsData = SimilarityWeights.getSharedModelData(config);

tos_data = modelDataFull["tos"]
models_out = ["EC-Earth3"]
SimilarityWeights.keepModelSubset(tos_data, findall(m -> !(m in models_out), Array(dims(tos_data, :model))))
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

for var in config.variables
    SimilarityWeights.keepModelSubset!(modelDataFull[var], shared_models);
    SimilarityWeights.keepModelSubset!(modelDataRef[var], shared_models);
end

weights = SimilarityWeights.overallWeights(
    modelDataRef, 
    obsData, 
    config.weight_contributions["performance"],
    config.weight_contributions["independence"], 
    config.weights_variables["performance"],
    config.weights_variables["independence"]   
);   

means = SimilarityWeights.getWeightedAverages(modelDataFull, weights);
