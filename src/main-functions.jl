import YAML
using DimensionalData


"""
    getSharedModelData(config::Config)

Return all data that is shared across variables in the reference period as well 
as in the full period.
"""
function getSharedModelData(config::Config)
    modelDataRef = loadDataFromConfig(config, config.name_ref_period, config.models_project_name);
    modelDataRef = getCommonModelsAcrossVars(modelDataRef);
    obsData = loadDataFromConfig(config, config.name_ref_period, config.obs_data_name);
    modelDataFull = loadDataFromConfig(config, config.experiment, config.models_project_name);
    modelDataFull = getCommonModelsAcrossVars(modelDataFull);
    
    # get weights just for those models in modelDataFull that are also in 
    # reference period 
    var = first(config.variables) # choose any variable
    diagnostic = first(config.diagnostics) # choose any diagnostic
    shared_models = intersect(
        modelDataFull[diagnostic][var].metadata["full_model_names"], 
        modelDataRef[diagnostic][var].metadata["full_model_names"]
    );
    for var in config.variables
        for diagnostic in config.diagnostics
            modelDataFull[diagnostic][var] = keepModelSubset(modelDataFull[diagnostic][var], shared_models);
            modelDataRef[diagnostic][var] = keepModelSubset(modelDataRef[diagnostic][var], shared_models);
        end
    end
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
    means = getWeightedAverages(modelDataFull["CLIM"], weights);

    models = weights.metadata["full_model_names"];
    model_key = getCMIPModelsKey(weights.metadata);
    @info "Nb included models (without ensemble members): " length(weights.metadata[model_key])
    foreach(m -> @info(m), models)
    saveWeights(weights, means, config.target_dir)
    return (weights=weights, avgs=means)
end
