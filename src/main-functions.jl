import YAML
using DimensionalData



function getSharedModelData(config::Config)
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
    return (modelDataFull, modelDataRef, obsData)
end

"""
    getWeightedAverages(
        modelDataAllVars::Dict{String, DimArray}, 
        weights::DimArray
    )
Compute average of 'modelDataAllVars', once weighted by 'weights' and once 
unweighted.
    
Note that the model dimension of 'weights' can be smaller than the model
dimension of 'modelDataAllVars' which may contain the predictions of all 
ensemble members. Here, these are averaged, s.t. for each model there is a 
just one prediction.
"""
function getWeightedAverages(modelDataAllVars::Dict{String, DimArray}, weights::DimArray)
    results = Dict{String, Dict{String, DimArray}}("weighted" => Dict(), "unweighted" => Dict());
    for var in keys(modelDataAllVars)
        data = modelDataAllVars[var];
        results["unweighted"][var] = computeWeightedAvg(data);
        results["weighted"][var] = computeWeightedAvg(data, weights);
    end
    return results
end

"""
    runWeights(path_config::String)

Compute weights and weighted/unweighted average of the data specified in the
config file located at 'path_config'.

"""
function runWeights(config::Config)
    modelDataFull, modelDataRef, obsData = getSharedModelData(config);
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
    saveWeights(weights, avgs, config.target_dir)
    return (weights=weights, avgs=means)
end
