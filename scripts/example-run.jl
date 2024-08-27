
using SimilarityWeights

path_config = "configs/example_historical_albedo.yml"
weights = SimilarityWeights.runWeights(path_config, true);


config = SimilarityWeights.validateConfig(path_config);
pathsDict = SimilarityWeights.buildPathsToVarData(config)
modelData = SimilarityWeights.loadPreprocData(pathsDict, [config.models_project_name]);
obsData = SimilarityWeights.loadPreprocData(pathsDict, [config.obs_data_name])

modelDataAllVars =  SimilarityWeights.getCommonModelsAcrossVars(modelData);