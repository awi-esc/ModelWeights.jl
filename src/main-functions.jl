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
    # TODO: this needs to be done differently
    shared_models = intersect(
        first(values(modelData1)).metadata["member_names"], 
        first(values(modelData2)).metadata["member_names"]
    );
    for id in keys(modelData1)
        modelData1[id] = subsetModelData(modelData1[id], shared_models)
    end 
    for id in keys(modelData2)       
        modelData2[id] = subsetModelData(modelData2[id], shared_models);
    end
end


"""
    computeWeights(model_data::Data, obs_data::Data, config::ConfigWeights)

Compute weight for each model in multi-model ensemble according to approach
from Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner,
Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming
from CMIP6 Projections When Weighting Models by Performance and
Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020):
995–1012. https://doi.org/10.5194/esd-11-995-2020.

# Arguments:
- `model_data`: Models for data for computing independence and performance weights.
- `obs_data`: Observational data for computing performance weights.
- `config`: parameters specifiying the relative contributions of each 
combination of variable and diagnostic.
"""
function computeWeights(
    model_data::Data, obs_data::Data, config_weights::ConfigWeights
)
    # only the observational data is used for which there is also model data
    obs_ids = filter(x -> x in model_data.ids, obs_data.ids)
    # and only those models are used for which there is observational data too
    model_ids = filter(x -> x in obs_ids, model_data.ids)
    # TODO: ids are not necessarily in same order, ordering should be defined
    if obs_ids != model_ids && !isempty(obs_ids)
        msg = "Model and observational data ids do not match! Obs: $obs_ids , Model: $model_ids"
        throw(ArgumentError(msg))
    end

    # make sure that weights (for diagnostic+variables) are normalized
    weights_perform = normalizeWeightsVariables(config_weights.performance)    
    weights_indep = normalizeWeightsVariables(config_weights.independence)

    distances_perform_all = computeDistancesAllDiagnostics(
        model_data, obs_data, collect(keys(config_weights.performance)), true
    )
    distances_indep_all = computeDistancesAllDiagnostics(
        model_data, obs_data, collect(keys(config_weights.independence)), false
    )
    Di = computeGeneralizedDistances(distances_perform_all, weights_perform, true)
    Sij = computeGeneralizedDistances(distances_indep_all, weights_indep, false)

    performances = performanceParts(Di, config_weights.sigma_performance)
    independences = independenceParts(Sij, config_weights.sigma_independence)
    weights = performances ./ independences;
    weights = weights ./ sum(weights);
    # TODO: add metadata
    #weights.metadata["name_ref_period"] = config_weights.ref_period  
    wP = performances ./ sum(performances)
    wI = independences ./ sum(independences)
    #w = wP./wI # just for sanity check

    climwip_weights =  ClimwipWeights(
        performance_distances = distances_perform_all,
        independence_distances = distances_indep_all, 
        Di = Di,
        Sij = Sij,
        wP = wP,
        wI = wI,
        w =  weights
        #overall = w./sum(w), # just for sanity check
    )
    logWeights(weights.metadata);
    if !isempty(config_weights.target_dir)
        saveWeights(climwip_weights, config_weights.target_dir)
    end
    return climwip_weights
end
