import SimilarityWeights as sw
using NCDatasets


#####  Example to compute weights for all data  #######
path_config = "configs/example_historical_albedo_weights.yml";
#path_config = "configs/example_historical_local.yml"

config_weights = sw.validateConfig(path_config);
weights = sw.getOverallWeights(confi_weights);

#### compute weighted averages ####
path_config = "configs/example_historical_albedo_graph.yml";
config = sw.validateConfig(path_config);
avgs = sw.getWeightedAverages(config, weights);

# TODO: add title weights based on which variables and diagnostics, add this to metadata when saving weights!
sw.plotWeights(weights)
sw.plotMeanData(config, avgs)

# convert to celsius
avgs["weighted"]["tas"] = sw.kelvinToCelsius(avgs["weighted"]["tas"]);
avgs["unweighted"]["tas"] = sw.kelvinToCelsius(avgs["unweighted"]["tas"]);
sw.plotMeanData(config, avgs)


### Look at performance and independence weights seperately ###
#w = sw.loadWeightsAsDimArray("/albedo/work/projects/p_forclima/britta/similarityweights-output/climwip-simplified/2024-09-24_14_32/weights.nc")
modelDataFull, modelDataRef = sw.loadModelData(config_weights);
obsData = sw.loadDataFromConfig(
    config_weights, 
    config_weights.name_obs_period, 
    config_weights.obs_data_name
);

wP = sw.getPerformanceWeights(modelDataRef["CLIM"], obsData["CLIM"], config_weights.weights_variables["performance"]);
# Di = sw.reduceGeneralizedDistancesVars(wP);
Di = sw.getPerformanceWeights(modelDataRef["CLIM"], obsData["CLIM"], config_weights.weights_variables["performance"], true);


wI = sw.getIndependenceWeights(modelDataRef["CLIM"], config_weights.weights_variables["independence"]);
# Sij = sw.reduceGeneralizedDistancesVars(wI);
Sij = sw.getIndependenceWeights(modelDataRef["CLIM"], config_weights.weights_variables["independence"], true);

figs_wP = sw.plotPerformanceWeights(wP)
figs_Di = sw.plotPerformanceWeights(wP, Di, false)
figs_Di = sw.plotPerformanceWeights(wP, nothing, false)

figs_wI = sw.plotIndependenceWeights(wI)
figs_Sij = sw.plotIndependenceWeights(Sij)


# build the actual weights from performance/independence weights
performances = sw.performanceParts(Di, config.weight_contributions["performance"]);
independences = sw.independenceParts(Sij, config.weight_contributions["independence"]);
sw.plotWeightContributions(independences, performances)

weights = performances ./ independences;
normalizedWeights = weights ./ sum(weights)


#########   Ensemble spread for some models with >1 ensemble member ########
data = modelDataRef["CLIM"]["tas"]
models_kept = ["EC-Earth3", "ACCESS-CM2"]
indices = findall(m -> m in models_kept, Array(dims(data, :model)))
shared_models = data.metadata["full_model_names"][indices]
data = sw.keepModelSubset(data, shared_models)
sw.plotEnsembleSpread(data, 7.5, 82.5)

###############################################################################
