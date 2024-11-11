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
ensemble members of that model. Instead of 'model', the returned DimArray has
dimension 'ensemble'.

# Arguments:
- `data`: a DimArray with at least dimension 'model'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
If true attribute ensemble_indices_map is set in metadata.
Set to false if vectors refer to different variables for instance. 
"""
function averageEnsembleVector(data::DimArray, updateMeta::Bool)
    data = renameModelDimsFromMemberToEnsemble(data, ["model"])
    
    grouped = groupby(data, :model=>identity);
    ensembles = collect(dims(grouped, :model))
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=:model), grouped)
    combined = cat(averages..., dims=(Dim{:model}(ensembles)));
    combined = set(combined, :model => :ensemble)

    if updateMeta
        meta = updateGroupedDataMetadata(data.metadata, grouped)
    else 
        meta = data.metadata
    end
    combined = rebuild(combined; metadata = meta);
    
    l = Lookups.Categorical(
        sort(ensembles);
        order=Lookups.ForwardOrdered()
    )
    combined = combined[ensemble=At(sort(ensembles))]
    combined = DimensionalData.Lookups.set(combined, ensemble=l)
    return combined
end


"""
    averageEnsembleMatrix(data::DimArray, updateMeta::Bool)

Compute the average weights across all members of each model ensemble for each
given variable.

# Arguments:
- `data`: DimArray with at least dimensions 'model1', 'model2'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
If true attribute ensemble_indices_map is set in metadata. 
Set to false if vectors refer to different variables for instance. 
"""
function averageEnsembleMatrix(data::DimArray, updateMeta::Bool)
    data = renameModelDimsFromMemberToEnsemble(data, ["model1", "model2"])
    
    ensembles = collect(unique(dims(data, :model1)))

    grouped = groupby(data, :model2=>identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=:model2), grouped)
    combined = cat(averages..., dims=(Dim{:model2}(ensembles)))

    grouped = groupby(combined, :model1=>identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=:model1), grouped)
    combined = cat(averages..., dims=(Dim{:model1}(ensembles)));
    
    for m in ensembles
        combined[model1=At(m), model2=At(m)] .= 0
    end
    if updateMeta
        meta = updateGroupedDataMetadata(data.metadata, grouped)
    else 
        meta = data.metadata
    end
    combined = rebuild(combined; metadata = meta);

    l = Lookups.Categorical(
        sort(ensembles);
        order=Lookups.ForwardOrdered()
    )
    combined = combined[model1=At(sort(ensembles)), model2=At(sort(ensembles))]
    combined = DimensionalData.Lookups.set(combined, model1=l, model2=l)
    combined = set(combined, :model1 => :ensemble1, :model2 => :ensemble2)
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
(i.e. ensemble, not ensemble members) pair has two dimensions, 'ensemble1' and 'ensemble2'
- `sigmaS`: free model parameter for impact of independence weights
"""
function independenceParts(generalizedDistances::DimArray, sigmaS::Number)
    # note: (+1 in eq. in paper is from when model is compared to itself since exp(0)=1)
    indep_parts = sum(exp.(-(generalizedDistances ./ sigmaS).^2), dims=:ensemble2);
    indep_parts = dropdims(indep_parts, dims=:ensemble2)
    indep_parts = set(indep_parts, :ensemble1 => :ensemble);
    return indep_parts
end


"""
    computeWeightedAvg(data::DimArray, w::Union{DimArray, Nothing}=nothing)

Compute the average values for each (lon,lat) grid point in 'data_var', weighted
by weights 'w'. Members of the same ensemble are averaged before. If no weight
vector is provided, unweighted average is computed.

# Arguments:
- `data`: has dimensions lon, lat, ensemble
- `w`: weights for ensembles (not individual members); has dimension 'ensemble'
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
        weights = DimArray(w, Dim{:ensemble}(Array(dims(data, :ensemble))))
    else
        # weights may have been computed wrt a different set of variables as we use here, 
        # so the list of models for which weights have been computed may be shorter 
        # than the models of the given data (for the same reference period).
        if sort(collect(dims(weights, :ensemble))) != sort(collect(dims(data, :ensemble)))
            @warn "Mismatch between models that weights were computed for and models in the data."
        end
        data = data[ensemble = Where(m -> m in dims(weights, :ensemble))]
        n_models_data = length(dims(data, :ensemble))
        n_weights = length(weights)
        if n_models_data != n_weights
            msg = "nb of models for observational and model predictions does not match: ";
            msg2 = "weights: " * string(n_weights) * " , data: " * string(n_models_data);
            throw(ArgumentError(msg * msg2))
        end
    end
    @assert isapprox(sum(weights), 1; atol=10^-4)
    for m in dims(data, :ensemble)
        data[ensemble = At(m)] = data[ensemble = At(m)] .* weights[ensemble = At(m)]
    end

    weighted_avg = dropdims(reduce(+, data, dims=:ensemble), dims=:ensemble)
    return weighted_avg
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


"""
    saveWeights(
        weights::DimVector,
        target_dir::String;
        target_fn::String="weights.nc"
    )

# Arguments:
- `weights`:
- `target_dir`:
- `target_fn`:
"""
function saveWeights(
    weights::ClimwipWeights,
    target_dir::String;
    target_fn::String=""
)
    if !isdir(target_dir)
        mkpath(target_dir)
    end
    if isempty(target_fn)
       path_to_target = joinpath(target_dir, join(["weights", getCurrentTime(), ".nc"], "_", ""))
    else
        path_to_target = joinpath(target_dir, target_fn);
        if isfile(path_to_target)
            msg1 = "File: " * path_to_target * " already exists!";
            path_to_target = joinpath(target_dir, join([getCurrentTime(), target_fn], "_"))
            msg2 = "Weights saved as: " * path_to_target
            @warn msg1 * msg2
        end
    end
    ds = NCDataset(path_to_target, "c")

    ensembles = map(x -> string(x), dims(weights.w, :ensemble))
    defDim(ds, "ensemble", length(ensembles))
    # Add a new variable to store the model names
    defVar(ds, "ensemble", Array(ensembles), ("ensemble",))
    
    function addNCDatasetVar!(ds, dimensions, name)
        defDim(ds, name, length(dimensions))
        defVar(ds, name, Array(dimensions), (name,))
        return nothing  
    end

    # add dimension variables
    addNCDatasetVar!(ds, dims(weights.performance_distances, :model), "model")
    addNCDatasetVar!(ds, dims(weights.performance_distances, :variable), "variable")
    addNCDatasetVar!(ds, dims(weights.performance_distances, :diagnostic), "diagnostic")

    # global attributes
    # TODO: check metadata, standard_name shouldnt be present!
    for (k, v) in weights.w.metadata
        if isa(v, String) # doesnt work for complex types like vectors or dicts!
            ds.attrib[k] = v
        end
    end

    # Add actual data weights
    for name in fieldnames(ClimwipWeights)
        if String(name) in ["w", "wP", "wI", "Di"]
            v = defVar(ds, String(name), Float64, ("ensemble",))
            data = getfield(weights, name)
            v[:] = Array(data)
        elseif String(name) == "independence_distances"
            v = defVar(ds, String(name), Float64, ("model", "model", "variable", "diagnostic"))
            v[:,:,:,:] =  Array(weights.independence_distances)
        elseif String(name) == "Sij"
            v = defVar(ds, String(name), Float64, ("ensemble", "ensemble"))
            v[:,:] = Array(weights.Sij)
        elseif String(name) == "performance_distances"
            v = defVar(ds, String(name), Float64, ("model", "variable", "diagnostic"))
            v[:,:,:] = Array(weights.performance_distances)
       end 
    end
    close(ds)

    @info "saved data to " path_to_target
end


function applyWeights(data::Data, weights::ClimwipWeights)
    
end