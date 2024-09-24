import SimilarityWeights as sw
using NCDatasets

path_config = "configs/example_historical_albedo.yml"
path_config = "configs/example_historical_local.yml"
path_config = "configs/climwip_simplified.yml"
config = sw.validateConfig(path_config);

######################## runWeights ###########################################
weights, avgs = sw.runWeights(config);

w = sw.loadWeightsAsDimArray("/albedo/work/projects/p_forclima/britta/similarityweights-output/climwip-simplified/2024-09-24_14_32/weights.nc")

sw.plotWeights(weights)
sw.plotMeanData(config, avgs)

####################### overallWeights step by step: ##############################
modelDataFull, modelDataRef, obsData = sw.getSharedModelData(config);
weights = sw.overallWeights(
    modelDataRef, 
    obsData, 
    config.weight_contributions["performance"],
    config.weight_contributions["independence"], 
    config.weights_variables["performance"],
    config.weights_variables["independence"]   
);
means = sw.getWeightedAverages(modelDataFull["CLIM"], weights);
wP = sw.generalizedDistancesPerformance(modelDataRef["CLIM"], obsData["CLIM"], config.weights_variables["performance"]);
Di = sw.reduceGeneralizedDistancesVars(wP);
wI = sw.generalizedDistancesIndependence(modelDataRef["CLIM"], config.weights_variables["independence"]);
Sij = sw.reduceGeneralizedDistancesVars(wI);

figs_wP = sw.plotPerformanceWeights(wP)
figs_Di = sw.plotPerformanceWeights(wP, Di, false)
figs_Di = sw.plotPerformanceWeights(wP, nothing, false)


figs_wI = sw.plotIndependenceWeights(wI)
figs_Sij = sw.plotIndependenceWeights(Sij)



performances = sw.performanceParts(Di, config.weight_contributions["performance"]);
independences = sw.independenceParts(Sij, config.weight_contributions["independence"]);
sw.plotWeightContributions(independences, performances)


weights = performances ./ independences;
normalizedWeights = weights ./ sum(weights)
###############################################################################

data = modelDataFull["tos"]
models_out = ["EC-Earth3"]
indices = findall(m -> !(m in models_out), Array(dims(data, :model)));
shared_models = data.metadata["full_model_names"][indices]
data = sw.keepModelSubset(data, shared_models)
sw.plotEnsembleSpread(data, 7.5, 82.5)

###############################################################################
