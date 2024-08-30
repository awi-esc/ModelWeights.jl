import YAML
using DimensionalData

"""
    runWeights(path_config::String)

Compute weights and weighted/unweighted average of the data specified in the
config file located at 'path_config'.

"""
function runWeights(config::Config)
    modelDataRef = loadDataFromConfig(config, "name_ref_period", "models_project_name");
    modelDataRef = getCommonModelsAcrossVars(modelDataRef);
    obsData = loadDataFromConfig(config, "name_ref_period", "obs_data_name");
    # TODO: make sure that there is observational data for all variables for 
    # the respective reference period
    weights = overallWeights(
        modelDataRef, 
        obsData, 
        config.weight_contributions["performance"],
        config.weight_contributions["independence"], 
        config.weights_variables["performance"],
        config.weights_variables["independence"]   
    );
        
    modelDataFull = loadDataFromConfig(config, "name_full_period", "models_project_name");
    modelDataFull = getCommonModelsAcrossVars(modelDataFull);
    # get Data just for those models for which also historical reference period is available
    models = getModelsInRefPeriod(modelDataFull, modelDataRef);
    means = getWeightedAverages(models, weights);
    return (weights=weights, avgs=means)
end
