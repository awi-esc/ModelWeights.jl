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
    areaWeightMatrix = repeat(areaWeights', length(DimensionalData.dims(sqDiff, :lon)), 1);  
    weightedValues = areaWeightMatrix .* sqDiff;

    #weightMatrix = repeat(areaWeights', length(DimensionalData.dims(m1, :lon)), 1);  
    areaWeightMatrix = ifelse.(mask .== 1, 0, areaWeightMatrix); 
    normalization = sum(areaWeightMatrix);

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
    nbModels = length(dims(data, :member));
    maskMissing = dropdims(any(ismissing, data, dims=:member), dims=:member);

    matrixS = zeros(nbModels, nbModels);    
    for (i, model_i) in enumerate(eachslice(data; dims=:member))
        model_i[maskMissing .== 1] .= 0;
        
        for (j, model_j) in enumerate(eachslice(data[:, :, i+1:end]; dims=:member))
            idx = j + i;
            model_j[maskMissing .== 1] .= 0;
            s_ij = areaWeightedRMSE(model_i, model_j, maskMissing);
            matrixS[i, idx] = s_ij
        end
    end
    # make symmetrical matrix
    symDistMatrix = matrixS .+ matrixS';
    dim = Array(dims(modelData, :member));
    return DimArray(symDistMatrix, (Dim{:member1}(dim), Dim{:member2}(dim)), metadata=modelData.metadata)
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
    member_names = Vector{String}();
    for (i, model_i) in enumerate(eachslice(models; dims=:member))
        name = dims(models, :member)[i]  # Access the name of the current model
        maskNbMissing = (ismissing.(observations) + ismissing.(model_i)) .> 0; # observations or model is missing (or both)
        maskedObs = deepcopy(observations);
        # maskedObs = dropdims(ifelse.(maskNbMissing .> 0, 0, maskedObs), dims=:source);
        maskedObs = ifelse.(maskNbMissing .> 0, 0, maskedObs)

        maskedModel = ifelse.(maskNbMissing .> 0, 0, model_i);
        mse = areaWeightedRMSE(maskedModel, maskedObs, maskNbMissing);

        push!(member_names, name);
        push!(distances, mse);
    end
    return DimArray(distances, (Dim{:member}(member_names)), metadata=models.metadata)
end


"""
    normalizeWeightsVariables(weights_dict::Dict{String, Number})

Normalize weights for combinations of variable and diagnostic such that they 
sum up to 1.

# Arguments:
- `weights_dict`: maps from VARIABLE_DIAGNOSTIC (e.g. tas_CLIM) to weights

# Return:
DimArray with dimensions 'variable' and 'diagnostic' containig the normalized 
weights.
"""
function normalizeWeightsVariables(weights_dict::Dict{String, Number})
    filter!(((k,v),) -> v != 0, weights_dict)
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
    summarizeEnsembleMembersVector(
        data::DimArray, updateMeta::Bool; fn::Function=Statistics.mean
)

For each model and variable (if several given), compute the mean across all
members of that model. Instead of 'member', the returned DimArray has
dimension 'model'.

# Arguments:
- `data`: a DimArray with at least dimension 'model'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
Set to false if vectors refer to different variables for instance. 
"""
function summarizeEnsembleMembersVector(
    data::DimArray, updateMeta::Bool; fn::Function=Statistics.mean
)
    data = setLookupsFromMemberToModel(data, ["member"])
    grouped = groupby(data, :model=>identity);
    models = String.(collect(dims(grouped, :model)))
    averages = map(entry -> mapslices(x -> fn(skipmissing(x)), entry, dims=:model), grouped)
    combined = cat(averages..., dims=(Dim{:model}(models)));

    meta = updateMeta ? updateGroupedDataMetadata(data.metadata, grouped) : data.metadata
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

Compute the average weights across all members of each model for each given variable.

# Arguments:
- `data`: DimArray with at least dimensions 'member1', 'member2'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
Set to false if vectors refer to different variables for instance. 
"""
function averageEnsembleMatrix(data::DimArray, updateMeta::Bool)
    data = setLookupsFromMemberToModel(data, ["member1", "member2"])
    models = String.(collect(unique(dims(data, :model1))))

    grouped = groupby(data, :model2=>identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=:model2), grouped)
    combined = cat(averages..., dims=(Dim{:model2}(models)))

    grouped = groupby(combined, :model1=>identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=:model1), grouped)
    combined = cat(averages..., dims=(Dim{:model1}(models)));
    
    for m in models
        combined[model1=At(m), model2=At(m)] .= 0
    end

    meta = updateMeta ? updateGroupedDataMetadata(data.metadata, grouped) : data.metadata
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
(i.e. not model members) pair has two dimensions, 'model1' and 'model2'
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
    computeWeightedAvg(data::DimArray; weights::Union{DimArray, Nothing}=nothing)

Compute the average values for each (lon,lat) grid point in 'data', weighted
by weights 'weights'. If no weight vector is provided, unweighted average is computed.

# Arguments:
- `data`: has dimensions lon, lat, model/member
- `weights`: weights for models or individual members; has dimension model/memnber
"""
function computeWeightedAvg(
    data::DimArray; weights::Union{DimArray, Nothing}=nothing
)
    data = deepcopy(data)
    dim_symbol = hasdim(data, :member) ? :member : :model
    models_data = collect(dims(data, dim_symbol))
    
    if isnothing(weights)
        weights = makeEqualWeights(data.metadata, dim_symbol)
        models_weights = collect(dims(weights, dim_symbol))
    else
        models_weights = collect(dims(weights, dim_symbol))
        if sort(models_data) != sort(models_weights)
            # weights may have been computed wrt a different set of variables as we use here, 
            # so the list of models for which weights have been computed may be shorter 
            # than the models of the given data. But for all given weights there must be data 
            # for now.
            data_no_weights = [model for model in models_data if !(model in models_weights)]           
            @warn "No weights were computed for these models which are therefore not considered in the weighted average:" data_no_weights
            if any(x -> !(x in models_data), models_weights)
                msg = "Data of models missing for which weights have been computed."
                throw(ArgumentError(msg))
            end
        end
        # subset data s.t only data that will be weighted is considered
        if dim_symbol == :member
            data = data[member = Where(m -> m in models_weights)]
        else
            data = data[model = Where(m -> m in models_weights)]
        end
    end
    @assert isapprox(sum(weights), 1; atol=10^-4)
    # there shouldnt be data for which there is no weight since it was filtered out above
    @assert Array(models_weights) == Array(dims(data, dim_symbol))
    
    # TODO: the following if-differentiation should not be necessary, find out how to 
    # index dimarray with variable instead of direct name
    if dim_symbol == :member
        for m in models_weights
            data[member = At(m)] = data[member = At(m)] .* weights[member = At(m)]
        end
    else
        for m in models_weights
            data[model = At(m)] = data[model = At(m)] .* weights[model = At(m)]
        end
    end
    weighted_avg = dropdims(reduce(+, data, dims=dim_symbol), dims=dim_symbol)
    return weighted_avg
end


"""
    makeEqualWeights(metadata::Dict, dimension::Symbol)

Equally distribute weight for each model over its members. Metadata is used
to access the individual model members the computed weights were based on.
    
# Arguments:
- `metadata`: metadata of the data for which weights were computed.
- `dimension`: level of data, e.g. 'member', 'model'
"""
function makeEqualWeights(metadata::Dict, dimension::Symbol)

    model_members = metadata["member_names"]
    n_models = length(model_members) # member_names is a vector of vectors! 
    if dimension == :member
        # make sure that the number of members per model is considered
        w = [repeat([1/n_models * (1/length(members))], length(members)) for members in model_members]
        w = vcat(w...)
        dimnames = vcat(model_members...)
    else
        w = [1/n_models for _ in range(1, n_models)]
        dimnames = unique(metadata["model_names"])
    end
    weights = DimArray(w, Dim{dimension}(dimnames), metadata = deepcopy(metadata))
    return weights
end


function logWeights(metadata_weights)
    @info "Nb included models (without members): " length(unique(metadata_weights["model_names"]))
    foreach(m -> @info(m), metadata_weights["member_names"])
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
    weights::ModelWeights,
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

    models = map(x -> string(x), dims(weights.w, :model))
    defDim(ds, "model", length(models))
    # Add a new variable to store the model names
    defVar(ds, "model", Array(models), ("model",))
    
    # add dimension variables
    function addNCDatasetVar!(ds, dimensions, name)
        defDim(ds, name, length(dimensions))
        defVar(ds, name, Array(dimensions), (name,))
        return nothing  
    end
    addNCDatasetVar!(ds, dims(weights.performance_distances, :member), "member")
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
    for name in fieldnames(ModelWeights)
        if String(name) in ["w", "wP", "wI", "Di"]
            v = defVar(ds, String(name), Float64, ("model",))
            data = getfield(weights, name)
            v[:] = Array(data)
        elseif String(name) == "independence_distances"
            v = defVar(ds, String(name), Float64, ("member", "member", "variable", "diagnostic"))
            v[:,:,:,:] =  Array(weights.independence_distances)
        elseif String(name) == "Sij"
            v = defVar(ds, String(name), Float64, ("model", "model"))
            v[:,:] = Array(weights.Sij)
        elseif String(name) == "performance_distances"
            v = defVar(ds, String(name), Float64, ("member", "variable", "diagnostic"))
            v[:,:,:] = Array(weights.performance_distances)
       end 
    end
    close(ds)
    @info "saved data to " path_to_target
end



""" 
    applyWeights(model_data::DimArray, weights::DimArray)

Compute the weighted average of model data 'data' with given weights 'weights'.
If the weights were computed for a subset of the models in 'data', they are normalized
and applied to the subset. Only weights per model (not members) are considered
for now, in the future, members should be considered too.

# Arguments:
- `model_data`: model predictions. If given for model members, the predictions of 
each model are considered the average value of all members of the respective model.
- `weights`: if given for each member of a model, these will be summed up to yield one
value per model.
"""
function applyWeights(model_data::DimArray, weights::DimArray)
    if hasdim(weights, :member)
        weights = summarizeEnsembleMembersVector(weights, true; fn=sum)
    end
    if hasdim(model_data, :member)
        # take average over model predictions of members of same model
        model_data = summarizeEnsembleMembersVector(model_data, true; fn=Statistics.mean)
    end
    models_weights = dims(weights, :model)
    models = dims(model_data, :model)
    shared_models = intersect(models_weights, models)

    if length(shared_models) == 0
        throw(ArgumentError("No models given for which weights had been computed!"))
    else 
        data_out = model_data[model = Where(x -> !(x in shared_models))]
        model_data = model_data[model = Where(x -> x in shared_models)]
        weights = weights[model = Where(x -> x in shared_models)]
        weights = weights ./ sum(weights)
        models_out = dims(data_out, :model)
        if length(models_out) > 0
            @warn "No weights had been computed for $models_out"
        end
        if length(shared_models) < length(models_weights)
            @warn "Weights were computed for a subset of the models of the given data. Weights were renormalized."
        end
    end
    return computeWeightedAvg(model_data; weights)
end