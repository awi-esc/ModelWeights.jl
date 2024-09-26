import SimilarityWeights as sw
using NCDatasets

begin
    path_config_weights = "configs/climwip_simplified_weights.yml";
    config_weights = sw.validateConfig(path_config_weights);

    # get overall weights
    weights = sw.computeWeights(config_weights);
    sw.plotWeights(weights)

    ##### Performance and Independence weights #####
    modelDataFull, modelDataRef = sw.loadModelData(config_weights);
    modelDataFull, modelDataRef = sw.getSharedModelData(modelDataFull, modelDataRef);
    obsData = sw.loadDataFromConfig(
        config_weights, 
        config_weights.name_ref_period, 
        config_weights.obs_data_name
    );

    wP = sw.generalizedDistancesPerformance(
        modelDataRef["CLIM"], 
        obsData["CLIM"], 
        config_weights.weights_variables["performance"]
    );
    Di = sw.reduceGeneralizedDistancesVars(wP);
    wI = sw.generalizedDistancesIndependence(
        modelDataRef["CLIM"], 
        config_weights.weights_variables["independence"]
    );
    Sij = sw.reduceGeneralizedDistancesVars(wI);

    figs_wP = sw.plotPerformanceWeights(wP)
    figs_Di = sw.plotPerformanceWeights(wP, Di, false)
    figs_Di = sw.plotPerformanceWeights(wP, nothing, false)

    figs_wI = sw.plotIndependenceWeights(wI)
    figs_Sij = sw.plotIndependenceWeights(Sij)
end

##### Temperature graph #####
    # plot of weighted averages of entire time period (1995-2014) 
    # avgs = sw.getWeightedAverages(config_weights, weights);
    # sw.plotMeanData(config_weights, avgs)
# TODO: double check results for temp. graph!
begin
    path_config_temp_graph = "configs/climwip_simplified_temperature_graph.yml";
    config_graph = sw.validateConfig(path_config_temp_graph);

    modelDataFull, modelDataRef = sw.loadModelData(config_graph);
    modelDataFull, modelDataRef= sw.getSharedModelData(modelDataFull, modelDataRef);

    tas_data = sw.averageEnsembleVector(modelDataFull["ANOM"]["tas"], true)
    sw.plotTempGraph(tas_data, weights, "1981-2010")
end

##### Temperature map #####
# TODO: add map plot here
begin
    path_config_temp_map = "configs/climwip_simplified_temperature_map.yml";
    config_map = sw.validateConfig(path_config_temp_map)
end
