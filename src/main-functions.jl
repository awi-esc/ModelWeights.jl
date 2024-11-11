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
    obs_keys = map(x -> (var = x.variable, stat = x.statistic), obs_data.ids)
    model_keys = map(x -> (var = x.variable, stat = x.statistic), model_data.ids)
    if sort(obs_keys) != sort(model_keys)
        msg = "Variable+Diagnostic combinations are not the same for observational and model data! Obs: $obs_keys , Model: $model_keys"
        throw(ArgumentError(msg))
    end
    # make sure that weights (for diagnostic+variables) are normalized
    weights_perform = getNormalizedWeightsVariables(config_weights.performance)
    weights_indep = getNormalizedWeightsVariables(config_weights.independence)

    # compute performance/independence distances for all models and ensemble members
    diagnostics = unique(map(x -> x.stat, obs_keys))
    Di = []
    Sij = []
    distances_perform_all = []
    distances_indep_all = []
    for diagnostic in diagnostics
        distances_perform = []
        distances_indep = []
        variables = unique(map(x -> x.var, model_keys))
        for var in variables
            k = var * "_" * diagnostic
            models_dict = filter(((key, val),) -> occursin(k, key), model_data.data)
            # check that only one dataset for combination of variable + diagnostic
            if length(models_dict) > 1
                @warn "more than one dataset for $var and $diagnostic in model data, first is taken!"
            end
            models = first(values(models_dict))
            
            obs_dict = filter(((key, val),) -> occursin(k, key), obs_data.data)
            if length(obs_dict) > 1
                @warn "more than one dataset for $var and $diagnostic in observational data, first is taken!"
            end
            observations = first(values(obs_dict))

            distsIndep = getModelDistances(models)
            push!(distances_indep, distsIndep)
            distsPerform = getModelDataDist(models, observations)
            push!(distances_perform, distsPerform)
        end
        distances_perform = cat(distances_perform..., dims = Dim{:variable}(collect(variables)));
        distances_indep = cat(distances_indep..., dims = Dim{:variable}(collect(variables)));
        push!(distances_perform_all, distances_perform)
        push!(distances_indep_all, distances_indep)
    end
    # put into DimArray
    distances_perform_all = cat(distances_perform_all..., dims = Dim{:diagnostic}(collect(diagnostics)));
    distances_indep_all = cat(distances_indep_all..., dims = Dim{:diagnostic}(collect(diagnostics)));
    
    # get normalizations for each diagnostic,variable combination: median across all models
    norm_perform = mapslices(Statistics.median, distances_perform_all, dims=:model)
    distances_perform = DimArray(
        distances_perform_all ./ norm_perform, 
        dims(distances_perform_all), 
        metadata = distances_perform_all.metadata
    )
    norm_indep = mapslices(Statistics.median, distances_indep_all, dims=(:model1, :model2))
    distances_indep = DimArray(
        distances_indep_all ./ norm_indep, 
        dims(distances_indep_all),
        metadata = distances_indep_all.metadata
    )

    #  take mean diagnostic per model (averaging ensemble members)
    distances_perform = averageEnsembleVector(distances_perform, false)
    distances_indep = averageEnsembleMatrix(distances_indep, false)
    
    # compute generalized distances Di, Sij (weighted average over computed distances)
    distances_perform = mapslices(x -> x .* weights_perform, 
        distances_perform, dims=(:variable, :diagnostic)
    )
    distances_indep = mapslices(x -> x .* weights_indep, 
        distances_indep, dims=(:variable, :diagnostic)
    )
    Di = dropdims(sum(distances_perform, dims=(:variable, :diagnostic)),
        dims=(:variable, :diagnostic)
    )
    Sij =  dropdims(sum(distances_indep, dims=(:variable, :diagnostic)),
        dims=(:variable, :diagnostic)
    )

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
