import YAML
using DimensionalData


"""
    keepSharedModelData(    
        modelDataFull::Dict{String, DimArray},
        modelDataRef::Dict{String, DimArray}
    )

Return the data of those models available in both datasets.
"""
function keepSharedModelData!(
    modelData1::Dict{String, DimArray},
    modelData2::Dict{String, DimArray}
)    
    shared_models = intersect(
        first(values(modelData1)).metadata["full_model_names"], 
        first(values(modelData2)).metadata["full_model_names"]
    );
    for id in keys(modelData1)
        modelData1[id] = getModelSubset(modelData1[id], shared_models)
    end 
    for id in keys(modelData2)       
        modelData2[id] = getModelSubset(modelData2[id], shared_models);
    end
end


"""
    getOverallWeights(model_data::Data, obs_data::Data, config::ConfigWeights)

Compute weight for each model in multi-model ensemble according to approach
from Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner,
Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming
from CMIP6 Projections When Weighting Models by Performance and
Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020):
995–1012. https://doi.org/10.5194/esd-11-995-2020.

# Arguments:
- `model_data`:
- `obs_data`:
- `config`:
"""
function getOverallWeights(
    model_data::Data, obs_data::Data, config::ConfigWeights
)::ClimwipWeights

    weights = computeWeights(model_data, obs_data, config);
    logWeights(weights.w.metadata);
    if !isempty(config.target_dir)
        saveWeights(weights, config.target_dir)
    end
    return weights
end


function getPerformanceWeights(
    modelData::Dict{String, DimArray}, 
    obsData::Dict{String, DimArray}, 
    weights_variables::Dict{String, Number}=Dict{String, Number}(),
    summarize_variables::Bool=false
)
    wP = computeDistances(
        modelData, "performance", weights_variables; obsData
    )
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
    wI = computeDistances(modelData, "independence", weights_variables);
    if summarize_variables
        wI = reduceGeneralizedDistancesVars(wI)
    end
    return wI
end


