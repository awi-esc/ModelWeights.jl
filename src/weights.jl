using NCDatasets
using DimensionalData
using Statistics
using LinearAlgebra

"""
    areaWeightedRMSE(m1::DimArray, m2::DimArray, mask::DimArray)

Compute the area weighted (cosine of latitudes in radians) root mean squared 
error between two matrices. 

# Arguments:
- `m1`: has dimensions 'lon', 'lat'
- `m2`: has dimensions 'lon', 'lat'
- `mask`: has values 0,1. Locations where mask is 1 are ignored, i.e. they get a weight of 0!

# Return:
- single value, area-weighted root mean squared error
"""
function areaWeightedRMSE(m1::DimArray, m2::DimArray, mask::DimArray)
    latitudes = dims(m1, :lat);
    areaWeights = cos.(deg2rad.(latitudes));
    areaWeights = DimArray(areaWeights, Dim{:lat}(Array(latitudes)));

    sqDiff = (m1 .- m2).^2;
    weightedValues = areaWeights' .* sqDiff;

    weightMatrix = repeat(areaWeights', length(DimensionalData.dims(m1, :lon)), 1);  
    weightMatrix = ifelse.(mask .== 1, 0, weightMatrix); 
    normalization = sum(weightMatrix);

    return sqrt(sum(skipmissing(weightedValues)./normalization));
end


"""
    getModelDistances(modelData::DimArray)

Compute the area weighted root mean squared error between model predictions for each pair of models. 

# Arguments:
- `modelData` has dimensions 'lon', 'lat', 'model' and contains the data for a single climate variable.

# Return:
- Symmetrical matrix (::DimArray) of size nxn where n is the number of models in 'modelData'. 
"""
function getModelDistances(modelData::DimArray)
    # Make sure to use a copy of the data, otherwise, it will be modified by applying the mask!!
    data = deepcopy(modelData);
    # only take values where none (!) of the models has infinite values!! (Not just the two that are compared to one another)
    nbModels = length(dims(data, :model));
    maskMissing = dropdims(any(ismissing, data, dims=:model), dims=:model);

    matrixS = zeros(nbModels, nbModels);    
    for (i, model_i) in enumerate(eachslice(data; dims=:model))
        model_i[maskMissing .== 1] .= 0;
        
        for (j, model_j) in enumerate(eachslice(data[:, :, i+1:end]; dims=:model))
            idx = j + i;
            model_j[maskMissing .== 1] .= 0;
            s_ij = areaWeightedRMSE(model_i, model_j, maskMissing);
            matrixS[i, idx] = s_ij
        end
    end
    # make symmetrical matrix
    symDistMatrix = matrixS .+ matrixS';
    dim = Array(dims(modelData, :model));
    return DimArray(symDistMatrix, (Dim{:model1}(dim), Dim{:model2}(dim)), metadata=modelData.metadata)
end


"""
    getModelDataDist(models::DimArray, observations::DimArray)

Compute the distance (area-weighted RMSE) between model predictions and observations. 

"""
function getModelDataDist(models::DimArray, observations::DimArray)      
    # Make sure to use a copy of the data, otherwise, it will be modified by applying the mask!!
    models = deepcopy(models);
    observations = deepcopy(observations);
    distances = [];
    model_names = Vector{String}();
    for (i, model_i) in enumerate(eachslice(models; dims=:model))
        model_name = dims(models, :model)[i]  # Access the name of the current model
        maskNbMissing = (ismissing.(observations) + ismissing.(model_i)) .> 0; # observations or model is missing (or both)
        maskedObs = deepcopy(observations);
        maskedObs = dropdims(ifelse.(maskNbMissing .> 0, 0, maskedObs), dims=:model);

        maskedModel = ifelse.(maskNbMissing .> 0, 0, model_i);
        mse = areaWeightedRMSE(maskedModel, maskedObs, maskNbMissing);

        push!(model_names, model_name);
        push!(distances, mse);
    end
    return DimArray(distances, (Dim{:model}(model_names)), metadata=models.metadata)
end


"""
    getNormalizedWeightsVariables(weightsVars::Dict{String, Number})

Normalize weights for each combination of variable and diagnostic, s.t. weights 
sum up to 1.

# Arguments:
- `weights_dict`: maps from VARIABLE_DIAGNOSTIC (e.g. tas_CLIM) to weights

# Return:
DimArray with dimensions 'variable' and 'diagnostic' containig the normalized 
weights.
"""
function getNormalizedWeightsVariables(weights_dict::Dict{String, Number})
    total = sum(values(weights_dict))
    normalized_weights = Dict{String, Number}()
    for key in keys(weights_dict)
        normalized_weights[key] = weights_dict[key] / total 
    end
    weights = map(x -> split(x, "_"), collect(keys(weights_dict)))
    variables = unique(map(first, weights))
    diagnostics = unique(map(x -> x[2], weights))

    w = DimArray(
        zeros(length(variables), length(diagnostics)),
        (Dim{:variable}(variables), Dim{:diagnostic}(diagnostics))
    )
    for (var, diag) in weights
        k = var * "_" * diag
        w[variable=At(var), diagnostic=At(diag)] = normalized_weights[k]
    end
    return w
end


""" 
    averageEnsembleVector(data::DimArray, updateMeta::Bool)

For each model and variable (if several given), compute the mean across all
ensemble members of that model.

# Arguments:
- `data`: a DimArray with dimension 'model' and possibly 'variable'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
If true attribute ensemble_indices_map is set in metadata.
Set to false if vectors refer to different variables for instance. 
"""
function averageEnsembleVector(data::DimArray, updateMeta::Bool)
    data = renameModelDimsFromMemberToEnsemble(data, ["model"])
    
    grouped = groupby(data, :model=>identity);
    models = collect(dims(grouped, :model))
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=:model), grouped)
    combined = cat(averages..., dims=(Dim{:model}(models)));

    if updateMeta
        meta = updateGroupedDataMetadata(data.metadata, grouped)
    else 
        meta = data.metadata
    end
    combined = rebuild(combined; metadata = meta);
    
    l = Lookups.Categorical(
        sort(models);
        order=Lookups.ForwardOrdered()
    )
    combined = combined[model=At(sort(models))]
    combined = DimensionalData.Lookups.set(combined, model=l)
    return combined
end


"""
    averageEnsembleMatrix(data::DimArray, updateMeta::Bool)

Compute the average weights across all members of each model ensemble for each
given variable.

# Arguments:
- `data`: DimArray with dimensions 'model1', 'model2' and possibly 'variable'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
If true attribute ensemble_indices_map is set in metadata. 
Set to false if vectors refer to different variables for instance. 
"""
function averageEnsembleMatrix(data::DimArray, updateMeta::Bool)
    data = renameModelDimsFromMemberToEnsemble(data, ["model1", "model2"])
    
    models = collect(unique(dims(data, :model1)))

    grouped = groupby(data, :model2=>identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=:model2), grouped)
    combined = cat(averages..., dims=(Dim{:model2}(models)))

    grouped = groupby(combined, :model1=>identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=:model1), grouped)
    combined = cat(averages..., dims=(Dim{:model1}(models)));
    
    for m in models
        combined[model1=At(m), model2=At(m)] .= 0
    end
    if updateMeta
        meta = updateGroupedDataMetadata(data.metadata, grouped)
    else 
        meta = data.metadata
    end
    combined = rebuild(combined; metadata = meta);

    l = Lookups.Categorical(
        sort(models);
        order=Lookups.ForwardOrdered()
    )
    combined = combined[model1=At(sort(models)), model2=At(sort(models))]
    combined = DimensionalData.Lookups.set(combined, model1=l, model2=l)
    return combined
end


"""
    performanceParts(generalizedDistances::DimArray, sigmaD::Number)

Compute the performance part (numerator) of the overall weight for each model.

# Arguments:
- `generalizedDistances`: contains generalized distances Di for each model
- `sigmaD`: free model parameter for impact of performance weights
"""
function performanceParts(generalizedDistances::DimArray, sigmaD::Number)
    return exp.(-(generalizedDistances ./ sigmaD).^2);
end


"""
    independenceParts(generalizedDistances::DimArray, sigmaS::Number)

Compute the independence part (denominator) of the overall weight for each model.

# Arguments:
- `generalizedDistances`: contains generalized distances S_{i,j} for each model
pair has two dimensions, 'model1' and 'model2'
- `sigmaS`: free model parameter for impact of independence weights
"""
function independenceParts(generalizedDistances::DimArray, sigmaS::Number)
    # note: (+1 in eq. in paper is from when model is compared to itself since exp(0)=1)
    indep_parts = sum(exp.(-(generalizedDistances ./ sigmaS).^2), dims=:model2);
    indep_parts = dropdims(indep_parts, dims=:model2)
    indep_parts = set(indep_parts, :model1 => :model);
    return indep_parts
end


"""
    computeWeightedAvg(data::DimArray, w::Union{DimArray, Nothing}=nothing)

Compute the average values for each (lon,lat) grid point in 'data_var', weighted
by weights 'w'. Members of the same ensemble are averaged before. If no weight
vector is provided, unweighted average is computed.

# Arguments:
- `data`: has dimensions lon, lat, model
- `w`: weights for ensembles (not individual members); has dimension 'model'
"""
function computeWeightedAvg(
    data::DimArray; weights::Union{DimArray, Nothing}=nothing
)
    #makeWeightPerEnsembleMember(weights);
    data = deepcopy(data)
    if isnothing(weights)
        # make sure that the number of ensemble members per model is considered
        ensemble_names = data.metadata["ensemble_names"]
        indices = []
        for model in unique(ensemble_names)
            push!(indices, findall(x -> x==model, ensemble_names))
        end
        n_ensembles = length(indices)
        w = []
        for positions in indices
            n_members = length(positions)
            for _ in range(1, n_members)
                push!(w, (1/n_ensembles) *  (1/n_members))
            end
        end
        weights = DimArray(w, Dim{:model}(Array(dims(data, :model))))
    else
        # weights may have been computed wrt a different set of variables as we use here, 
        # so the list of models for which weights have been computed may be shorter 
        # than the models of the given data (for the same reference period).
        if sort(collect(dims(weights, :model))) != sort(collect(dims(data, :model)))
            @warn "Mismatch between models that weights were computed for and models in the data."
        end
        data = data[model = Where(m -> m in dims(weights, :model))]
        n_models_data = length(dims(data, :model))
        n_weights = length(weights)
        if n_models_data != n_weights
            msg = "nb of models for observational and model predictions does not match: ";
            msg2 = "weights: " * string(n_weights) * " , data: " * string(n_models_data);
            throw(ArgumentError(msg * msg2))
        end
    end
    @assert isapprox(sum(weights), 1; atol=10^-4)
    for m in dims(data, :model)
        data[model = At(m)] = data[model = At(m)] .* weights[model = At(m)]
    end

    weighted_avg = dropdims(reduce(+, data, dims=:model), dims=:model)
    return weighted_avg
end


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

    return ClimwipWeights(
        performance_distances = distances_perform_all,
        independence_distances = distances_indep_all, 
        Di = Di,
        Sij = Sij,
        wP = wP,
        wI = wI,
        w =  weights
        #overall = w./sum(w), # just for sanity check
    )
end

"""
    makeWeightPerEnsembleMember(weights::DimArray)


"""
function makeWeightPerEnsembleMember(weights::DimArray)
    full_model_names = weights.metadata["full_model_names"]
    models = weights.metadata[getCMIPModelsKey(weights.metadata)]
    nb_ensemble_members = [];
    for model in models
        n = count(member -> startswith(member, model), full_model_names)
        push!(nb_ensemble_members, n)
    end
    weights_all_members = []
    for (idx, w) in enumerate(weights)
        n_members = nb_ensemble_members[idx]
        for _ in range(1, n_members)
            push!(weights_all_members, w/nb_ensemble_members[idx])
        end
    end
    return DimArray(weights_all_members, Dim{:model}(full_model_names), metadata = weights.metadata)
end


function logWeights(metadata_weights)
    models = metadata_weights["full_model_names"];
    @info "Nb included models (without ensemble members): " length(metadata_weights["ensemble_names"])
    foreach(m -> @info(m), models)
end


function applyWeights(data::Data, weights::ClimwipWeights)
    
end