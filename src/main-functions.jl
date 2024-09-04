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
    # TODO: make sure that there is observational data for all variables for 
    # the respective reference period, for now this is just assumed but in 
    # some rare cases it may be wrong
    obsData = loadDataFromConfig(config, "name_ref_period", "obs_data_name");
    modelDataFull = loadDataFromConfig(config, "name_full_period", "models_project_name");
    modelDataFull = getCommonModelsAcrossVars(modelDataFull);
    
    # get weights just for those models in modelDataFull that are also in 
    # reference period 
    shared_models = intersect(
        first(values(modelDataFull)).metadata["full_model_names"], 
        first(values(modelDataRef)).metadata["full_model_names"]
    );
    keepModelSubset!(modelDataFull, shared_models);
    keepModelSubset!(modelDataRef, shared_models);
    
    weights = overallWeights(
        modelDataRef, 
        obsData, 
        config.weight_contributions["performance"],
        config.weight_contributions["independence"], 
        config.weights_variables["performance"],
        config.weights_variables["independence"]   
    );   
    means = getWeightedAverages(modelDataFull, weights);

    models = weights.metadata["full_model_names"];
    @info "Nb included models (without ensemble members): " length(weights.metadata["source_id"])
    foreach(m -> @info(m), models)
    return (weights=weights, avgs=means)
end
