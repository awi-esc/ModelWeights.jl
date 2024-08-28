
using SimilarityWeights

path_config = "configs/example_historical_albedo.yml"
config = SimilarityWeights.validateConfig(path_config);

weights, avgs = SimilarityWeights.runWeights(config, true);

SimilarityWeights.plotMeanData(config, avgs)
