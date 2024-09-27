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
    ref_period_weights = config.name_ref_period
    # if isempty(ref_period_weights)
    #     config.name_ref_period = weights.metadata["name_ref_period"]
    # else 
    if !isempty(ref_period_weights) && (ref_period_weights != weights.metadata["name_ref_period"])
        msg = "weights were computed for period: " * weights.metadata["name_ref_period"]; 
        msg2 = ", but in config it is set to" * ref_period_weights;
        throw(ArgumentError(msg * msg2))
    end
    # end
    modelDataFull, _ = loadModelData(config);
    
    results = Dict{String, Dict{String, DimArray}}("weighted" => Dict(), "unweighted" => Dict());
    # TODO: maybe dont hard code CLIM here
    climatologies = modelDataFull["CLIM"]
    for var in keys(climatologies)
        data = climatologies[var];
        results["unweighted"][var] = computeWeightedAvg(data);
        results["weighted"][var] = computeWeightedAvg(data, weights);
    end
    return results
end


"""
    getOverallWeights(config::Config)

Compute weight for each model in multi-model ensemble according to approach
from Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner,
Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming
from CMIP6 Projections When Weighting Models by Performance and
Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020):
995–1012. https://doi.org/10.5194/esd-11-995-2020.
"""
function getOverallWeights(config::Config)
    _, modelDataRef = loadModelData(config);
    obsData = loadDataFromConfig(
        config, 
        config.name_obs_period, 
        config.obs_data_name
    );
    weights = computeWeights(
        modelDataRef, 
        obsData, 
        config.weights_variables["performance"],
        config.weights_variables["independence"],
        config.weight_contributions["performance"],
        config.weight_contributions["independence"]
    )
    weights.metadata["name_ref_period"] = config.name_ref_period
    logWeights(weights.metadata);
    saveWeights(weights, config.target_dir)
    return weights
end


function getPerformanceWeights(
    modelData::Dict{String, DimArray}, 
    obsData::Dict{String, DimArray}, 
    weights_variables::Dict{String, Number}=Dict{String, Number}(),
    summarize_variables::Bool=false
)
    wP = generalizedDistancesPerformance(modelData, obsData, weights_variables)
    if summarize_variables
        wP = reduceGeneralizedDistancesVars(wP)
    end
    return wP
end


function getIndependenceWeights(
    modelData::Dict{String, DimArray}, 
    weights_variables::Dict{String, Number}=Dict{String, Number}(), 
    summarize_variables::Bool=false
)
    wI = generalizedDistancesIndependence(modelData, weights_variables);
    if summarize_variables
        wI = reduceGeneralizedDistancesVars(wI)
    end
    return wI
end


