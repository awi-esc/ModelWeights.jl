import YAML
using DimensionalData


"""
    getSharedModelData(config::Config)

Return all data that is shared across variables in the reference period as well 
as in the full period.
"""
function getSharedModelData(
    modelDataFull::Union{Nothing, Dict{String, Dict{String, DimArray}}}, 
    modelDataRef::Union{Nothing, Dict{String, Dict{String, DimArray}}}
)    
    if !any(isnothing.([modelDataFull, modelDataRef]))
        # get weights just for those models in modelDataFull that are also in 
        # reference period
        diagnostics = keys(modelDataFull)
        diagnostic = first(diagnostics) # choose any diagnostic
        variables = keys(modelDataFull[diagnostic])
        var = first(variables) # choose any variable
        shared_models = intersect(
            modelDataFull[diagnostic][var].metadata["full_model_names"], 
            modelDataRef[diagnostic][var].metadata["full_model_names"]
        );
        for var in variables
            for diagnostic in diagnostics
                modelDataFull[diagnostic][var] = keepModelSubset(modelDataFull[diagnostic][var], shared_models);
                modelDataRef[diagnostic][var] = keepModelSubset(modelDataRef[diagnostic][var], shared_models);
            end
        end
    end
    return (modelDataFull, modelDataRef)
end

"""
    getWeightedAverages(
        config::Config,
        weights::DimArray
    )
Compute average of 'modelDataAllVars', once weighted by 'weights' and once 
unweighted.
    
Note that the model dimension of 'weights' can be smaller than the model
dimension of 'modelDataAllVars' which may contain the predictions of all 
ensemble members. Here, these are averaged, s.t. for each model there is a 
just one prediction.
"""
function getWeightedAverages(config::Config, weights::DimArray)
    modelDataFull, modelDataRef = loadModelData(config);
    modelDataFull, modelDataRef = getSharedModelData(modelDataFull, modelDataRef);
    
    results = Dict{String, Dict{String, DimArray}}("weighted" => Dict(), "unweighted" => Dict());
    # TODO: maybe dont hard code CLIM here
    for var in keys(modelDataFull["CLIM"])
        data = modelDataAllVars[var];
        results["unweighted"][var] = computeWeightedAvg(data);
        results["weighted"][var] = computeWeightedAvg(data, weights);
    end
    return results
end

"""
    computeWeights(path_config::String)

Compute weights for the data specified in config file. According to climwip method 
from Brunner et al.
"""
function computeWeights(config::Config)
    modelDataFull, modelDataRef = loadModelData(config);
    modelDataFull, modelDataRef = getSharedModelData(modelDataFull, modelDataRef);
    obsData = loadDataFromConfig(
        config, 
        config.name_obs_period, 
        config.obs_data_name
    );
    weights = overallWeights(
        modelDataRef, 
        obsData, 
        config.weight_contributions["performance"],
        config.weight_contributions["independence"], 
        config.weights_variables["performance"],
        config.weights_variables["independence"]   
    );   
    models = weights.metadata["full_model_names"];
    model_key = getCMIPModelsKey(weights.metadata);
    @info "Nb included models (without ensemble members): " length(weights.metadata[model_key])
    foreach(m -> @info(m), models)
    saveWeights(weights, config.target_dir)
    return weights
end
