import YAML
using DimensionalData


"""
    keepSharedModelData(    
        modelDataFull::Dict{String, DimArray},
        modelDataRef::Dict{String, DimArray}
    )

Return all data that is shared across variables in the reference period as well 
as in the period of interested (refered to as full period).
"""
function keepSharedModelData!(
    modelDataFull::Dict{String, DimArray},
    modelDataRef::Dict{String, DimArray}
)    
    shared_models = intersect(
        first(values(modelDataFull)).metadata["full_model_names"], 
        first(values(modelDataRef)).metadata["full_model_names"]
    );
    for id in keys(modelDataFull)
        modelDataFull[id] = getModelSubset(modelDataFull[id], shared_models)
    end 
    for id in keys(modelDataRef)       
        modelDataRef[id] = getModelSubset(modelDataRef[id], shared_models);
    end
end


"""
    getOverallWeights(data::Data, config::ConfigWeights)

Compute weight for each model in multi-model ensemble according to approach
from Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner,
Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming
from CMIP6 Projections When Weighting Models by Performance and
Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020):
995–1012. https://doi.org/10.5194/esd-11-995-2020.

# Arguments:
- `data`:
- `config`:
"""
function getOverallWeights(data::Data, config::ConfigWeights)::ClimwipWeights
    weights = computeWeights(data, config);
    logWeights(weights.overall.metadata);
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
    wP = generalizedDistances(
        modelData, "performance"; weights_variables, obsData
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
    wI = generalizedDistances(modelData, "independence"; weights_variables);
    if summarize_variables
        wI = reduceGeneralizedDistancesVars(wI)
    end
    return wI
end


