
using SimilarityWeights

path_config = "configs/example_historical_albedo.yml"
#path_config = "configs/example_historical_local.yml"

config = SimilarityWeights.validateConfig(path_config);
weights, avgs = SimilarityWeights.runWeights(config);

SimilarityWeights.plotMeanData(config, avgs)



# pathsDictRef = SimilarityWeights.buildPathsToVarData(config, config.name_ref_period);
# modelDataRef = SimilarityWeights.loadPreprocData(pathsDictRef, [config.models_project_name]);
# obsData = SimilarityWeights.loadPreprocData(pathsDictRef, [config.obs_data_name]);
# modelDataAllVarsRef =  SimilarityWeights.getCommonModelsAcrossVars(modelDataRef);

# wP = SimilarityWeights.getPerformanceWeights(modelDataAllVarsRef, obsData, config.weights_variables["performance"]);
# wI = SimilarityWeights.getIndependenceWeights(modelDataAllVarsRef, config.weights_variables["independence"]);
# sigmas = config.weight_contributions
# weights = SimilarityWeights.combineWeights(wP, wI, sigmas["performance"], sigmas["independence"]);

# # use weights to compute weighted averages
# pathsDictFull = SimilarityWeights.buildPathsToVarData(config, config.name_full_period);
# modelDataFull = SimilarityWeights.loadPreprocData(pathsDictFull, [config.models_project_name]);
# modelDataAllVarsFull =  SimilarityWeights.getCommonModelsAcrossVars(modelDataFull);


