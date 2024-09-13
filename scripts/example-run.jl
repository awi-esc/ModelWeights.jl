
import SimilarityWeights as sw

path_config = "configs/example_historical_albedo.yml"
path_config = "configs/example_historical_local.yml"
config = sw.validateConfig(path_config);

######################## runWeights ###########################################
weights, avgs = sw.runWeights(config);
sw.plotWeights(weights)
sw.plotMeanData(config, avgs)

modelDataFull, modelDataRef, obsData = sw.getSharedModelData(config);
weights = sw.overallWeights(
    modelDataRef, 
    obsData, 
    config.weight_contributions["performance"],
    config.weight_contributions["independence"], 
    config.weights_variables["performance"],
    config.weights_variables["independence"]   
);      
means = sw.getWeightedAverages(modelDataFull, weights);
####################### overallWeights step by step: ##############################
weightsVarsPerform = Dict{String, Number}("tas" => 1, "pr" => 2, "tos" => 1); 
weightsVarsIndep = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "tos" => 0.25); 

wP = sw.generalizedDistancesPerformance(modelDataRef, obsData, weightsVarsPerform);
Di = sw.overallGeneralizedDistances(wP);
wI = sw.generalizedDistancesIndependence(modelDataRef, weightsVarsIndep);
Sij = sw.overallGeneralizedDistances(wI);

performances = sw.performanceParts(Di, 0.5);
independences = sw.independenceParts(Sij, 0.5);
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
