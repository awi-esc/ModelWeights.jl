
using SimilarityWeights

path_config = "configs/example_historical_albedo.yml"
#path_config = "configs/example_historical_local.yml"

config = SimilarityWeights.validateConfig(path_config);
weights, avgs = SimilarityWeights.runWeights(config);

SimilarityWeights.plotMeanData(config, avgs)
