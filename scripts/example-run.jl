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

"""
    plotMeanData(config::Config, means::Dict{String, Dict{String, DimArray}})

# Arguments:
- `config`: path to configuration file
- `means`: maps from 'weighted'/'unweighted' to climate variables to mean data.
"""
function plotMeanData(config::Config, means::Dict{String, Dict{String, DimArray}})
    for avg_type in keys(means)
        data = means[avg_type]
        for var in keys(data)
            var_data = data[var]
            title = join(
                [avg_type, var_data.metadata["long_name"], "in", 
                 var_data.metadata["units"], "\n experiment:",
                 var_data.metadata["experiment_id"]
                ], " ", " "
            );
            f = plotMeansOnMap(var_data,  title);
            sw.savePlot(f, config.target_dir, join([avg_type, var_data.metadata["standard_name"], "png"], "_", "."))
        end
    end
end
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
