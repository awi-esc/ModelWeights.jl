using DimensionalData
using NCDatasets
using Statistics
using LinearAlgebra
using JLD2
using Setfield
using Distributions


"""
    areaWeightedRMSE(m1::AbstractArray, m2::AbstractArray, mask::AbstractArray)

Compute the area weighted (cosine of latitudes in radians) root mean squared 
error between two matrices. 

# Arguments:
- `m1`: has dimensions 'lon', 'lat'
- `m2`: has dimensions 'lon', 'lat'
- `mask`: has values 0,1. Locations where mask is 1 are ignored, i.e. they get a weight of 0!

# Return:
- single value, area-weighted root mean squared error
"""
function areaWeightedRMSE(m1::AbstractArray, m2::AbstractArray, mask::AbstractArray)
    # TODO: update function, just use single input data Array and use new function 
    # to compute the weights
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
    getModelDistances(modelData::AbstractArray)

Compute the area weighted root mean squared error between model predictions for each pair of models. 

# Arguments:
- `modelData` has dimensions 'lon', 'lat', 'model' and contains the data for a single climate variable.

# Return:
- Symmetrical matrix (::YAXArray) of size nxn where n is the number of models in 'modelData'. 
"""
function getModelDistances(modelData::YAXArray)
    # Make sure to use a copy of the data, otherwise, it will be modified by applying the mask!!
    data = deepcopy(modelData);
    # only take values where none (!) of the models has infinite values!! (Not just the two that are compared to one another)
    nbModels = length(dims(data, :member));
    maskMissing = dropdims(any(ismissing, data, dims=:member), dims=:member);

    matrixS = zeros(nbModels, nbModels);    
    for (i, model_i) in enumerate(eachslice(data[:,:,1:end-1]; dims=:member))
        model_i_mat = Array(model_i)
        model_i_mat[maskMissing .== 1] .= 0;
        for (j, model_j) in enumerate(eachslice(data[:, :, i+1:end]; dims=:member))
            model_j_mat = Array(model_j)
            idx = j + i;
            model_j_mat[maskMissing .== 1] .= 0;
            s_ij = areaWeightedRMSE(model_i, model_j, maskMissing);
            matrixS[i, idx] = s_ij
        end
    end
    # make symmetrical matrix
    symDistMatrix = matrixS .+ matrixS';
    dim = Array(dims(modelData, :member));
    return YAXArray((Dim{:member1}(dim), Dim{:member2}(dim)), symDistMatrix, deepcopy(modelData.properties))
end


"""
    getModelDataDist(models::AbstractArray, observations::AbstractArray)

Compute the distance (area-weighted RMSE) between model predictions and observations. 

# Arguments:
- `models`:
- `observations`:
"""
function getModelDataDist(models::YAXArray, observations::AbstractArray)      
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
    return YAXArray((Dim{:member}(member_names),), distances, models.properties)
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
    averageEnsembleMembersMatrix(data::AbstractArray, updateMeta::Bool)

Compute the average across all members of each model for each given variable 
for model to model data, e.g. distances between model pairs.

# Arguments:
- `data`: with at least dimensions 'member1', 'member2'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
Set to false if vectors refer to different variables for instance. 
"""
function averageEnsembleMembersMatrix(data::YAXArray, updateMeta::Bool)
    data = setLookupsFromMemberToModel(data, ["member1", "member2"])
    models = String.(collect(unique(dims(data, :model1))))

    grouped = groupby(data, :model2=>identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=(:model2,)), grouped)
    combined = cat(averages..., dims=(Dim{:model2}(models)))

    grouped = groupby(combined, :model1=>identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims=(:model1,)), grouped)
    combined = cat(averages..., dims=(Dim{:model1}(models)));
    
    for m in models
        combined[model1=At(m), model2=At(m)] .= 0
    end

    meta = updateMeta ? updateGroupedDataMetadata(data.properties, grouped) : data.properties
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
    performanceParts(generalizedDistances::AbstractArray, sigmaD::Number)

Compute the performance part (numerator) of the overall weight for each model.

# Arguments:
- `generalizedDistances`: contains generalized distances Di for each model
- `sigmaD`: free model parameter for impact of performance weights
"""
function performanceParts(generalizedDistances::AbstractArray, sigmaD::Number)
    return exp.(-(generalizedDistances ./ sigmaD).^2);
end


"""
    independenceParts(generalizedDistances::AbstractArray, sigmaS::Number)

Compute the independence part (denominator) of the overall weight for each model.

# Arguments:
- `generalizedDistances`: contains generalized distances S_{i,j} for each model 
(i.e. not model members) pair has two dimensions, 'model1' and 'model2'
- `sigmaS`: free model parameter for impact of independence weights
"""
function independenceParts(generalizedDistances::AbstractArray, sigmaS::Number)
    # note: (+1 in eq. in paper is from when model is compared to itself since exp(0)=1)
    indep_parts = sum(exp.(-(generalizedDistances ./ sigmaS).^2), dims=:model2);
    indep_parts = dropdims(indep_parts, dims=:model2)
    indep_parts = DimensionalData.set(indep_parts, :model1 => :model);
    return indep_parts
end


"""
    computeWeightedAvg(
        data::YAXArray; 
        weights::Union{DimArray, Nothing}=nothing,
        use_members_equal_weights::Bool=true
    )

Compute the average values for each (lon,lat) grid point in 'data', weighted
by weights 'weights'. If no weight vector is provided, unweighted average is computed.

# Arguments:
- `data`: must have dimension 'model' or 'member'
- `weights`: weights for models or individual members
- `use_members_equal_weights`:  if `weights` is nothing, all models receive 
equal weight. If `use_members_equal_weights` is true, the number of members 
per model is taken into account, s.t. each model receives equal weight, which
is distributed among the respective members.
"""
function computeWeightedAvg(
    data::YAXArray; 
    weights::Union{DimArray, Nothing}=nothing, 
    use_members_equal_weights::Bool=true
)
    da = DimArray(Array{Union{Missing, Float64}}(undef, size(data)), dims(data))
    da .= Array(data)
    dim_symbol = hasdim(data, :member) ? :member : :model
    models_data = collect(dims(data, dim_symbol))
    
    if isnothing(weights)
        weights = makeEqualWeights(
            data.properties, dim_symbol, use_members_equal_weights
        )
        models_weights = collect(dims(weights, dim_symbol))
    else
        if !hasdim(weights, dim_symbol)
            msg = "level of weights and data must align (both :member or both :model)!"
            throw(ArgumentError(msg))
        end
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
            da = da[member = Where(m -> m in models_weights)]
        else
            da = da[model = Where(m -> m in models_weights)]
        end
    end
    @assert isapprox(sum(weights), 1; atol=10^-4)
    # there shouldnt be data for which there is no weight since it was filtered out above
    @assert sort(Array(models_weights)) == sort(Array(dims(data, dim_symbol)))
        
    # readjust weights for missing values, if one model has a missing value, 
    # it is ignored in the weighted average for that particular lon,lat-position.
    not_missing_vals = mapslices(x -> (ismissing.(x).==0), da; dims=(dim_symbol,))
    w_temp = mapslices(x -> (x .* weights) ./ sum(x .* weights), not_missing_vals, dims=(dim_symbol,))
    w_temp = replace(w_temp, NaN => missing)

    for m in models_weights
        value = getAtModel(da, dim_symbol, m) .* getAtModel(w_temp, dim_symbol, m)
        putAtModel!(da, dim_symbol, m, value)
    end        
    weighted_avg = mapslices(x -> sum(skipmissing(x)), da, dims=(dim_symbol,))
    weighted_avg = dropdims(weighted_avg, dims=dim_symbol)
    
    # set to missing when value was missing for ALL models
    n_models = length(dims(da, dim_symbol))
    all_missing = dropdims(
        mapslices(x -> sum(ismissing.(x)) == n_models, da, dims=(dim_symbol,)),
        dims = dim_symbol
    )
    if !any(all_missing)
        result = weighted_avg
    else
        result = Array{Union{Missing, Float64}}(undef, size(weighted_avg))
        s = repeat([:], length(size(weighted_avg)))
        result[s...] = Array(weighted_avg)
        result[all_missing] .= missing
    end
    return YAXArray(dims(weighted_avg), result, deepcopy(data.properties))
end


"""
    makeEqualWeights(metadata::Dict, dimension::Symbol; use_members::Bool=true)

Create a weight vector, with equal weight for each model. Distribute weight across
model members if dimension=:member. Metadata is used to access the individual 
model members the computed weights were based on.

# Arguments:
- `metadata`: metadata of the data for which weights were computed.
- `dimension`: level of data, e.g. 'member', 'model'
"""
function makeEqualWeights(
    metadata::Dict, dimension::Symbol, use_members::Bool=true
)
    model_members = metadata["member_names"]
    n_models = length(unique(metadata["model_names"]))
    if dimension == :member
        dimnames = model_members
        if use_members
            # make sure that the number of members per model is considered
            w = [repeat([1/n_models * (1/length(members))], length(members)) 
                for members in model_members]
            w = vcat(w...)
        else
            n_members = length(dimnames)
            w = repeat([1/n_members], n_members)
        end
    elseif dimension == :model
        w = [1/n_models for _ in range(1, n_models)]
        dimnames = unique(metadata["model_names"])
    else
        throw(ArgumentError("dimension must be one of: :model, :member"))
    end

    l = Lookups.Categorical(
        sort(dimnames);
        order=Lookups.ForwardOrdered()
    )
    weights = DimArray(w, Dim{dimension}(dimnames), metadata = deepcopy(metadata))
    if dimension == :member
        weights = weights[member = At(sort(dimnames))]
        weights = DimensionalData.Lookups.set(weights, member=l)
    else
        weights = weights[model = At(sort(dimnames))]
        weights = DimensionalData.Lookups.set(weights, model=l)
    end
    return weights
end


"""
    distributeWeightsAcrossMembers(weights::YAXArray)

Equally distribute weight for each model over its members. Metadata is used
to access the individual model members the computed weights were based on.

# Arguments:
- `weights`:
"""
function distributeWeightsAcrossMembers(weights::YAXArray)
    models = Array(dims(weights, :model))
    members = weights.properties["member_names"]
    w_members = []
    for (i, m) in enumerate(models)
        n_members = sum(weights.properties["model_names"] .== m)
        w = only(weights[model = At(m)])
        for _ in range(1, n_members) 
            push!(w_members, w/n_members)
        end
    end
    return YAXArray((Dim{:member}(members),), w_members, weights.properties)
end


function validateConfigTargetPath(target_path::String)
    target_dir = dirname(target_path)
    if !isdir(target_dir)
        mkpath(target_dir)
    end
    if isfile(target_path)
        msg = "$target_path already exisits, will use "
        target_path = joinpath(
            target_dir, join([getCurrentTime(), basename(target_path)], "_")
        )
        @info msg * "$target_path instead."
    end
    return target_path
end



"""
    saveWeightsAsNCFile(weights::Weights; target_path::String)

# Arguments:
- `weights`: Weights object to be saved.
- `target_path`: Path to where weights shall be stored.
"""
function saveWeightsAsNCFile(weights::Weights, target_path::String)
    target_path = validateConfigTargetPath(target_path)
    ds = NCDataset(target_path, "c")
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

    # Add actual data weights
    for name in fieldnames(Weights)
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
        elseif String(name) == "config"
            # configuration is added as global attributes:
            for config_name in fieldnames(ConfigWeights)
                config_val = getfield(weights.config, config_name)
                config_key = String(config_name)
                if config_key in ["performance", "independence"]
                    for (k,v) in config_val
                        ds.attrib["w_" * config_key * "_" * k] = v
                    end
                elseif config_key != "target_path"
                    # if target path was updated, use updated version here, 
                    # not the one that had been saved in weights!
                    ds.attrib[config_key] = config_val
                end
            end
       end 
    end
    close(ds)
    @info "saved weights to $(target_path)"
    return target_path
end


function saveWeightsAsJuliaObj(weights::Weights, target_path::String)
    target_path = validateConfigTargetPath(target_path)
    config = weights.config
    config = @set config.target_path = target_path
    weights = @set weights.config = config
    jldsave(target_path; weights=weights)
    @info "saved weights to: $(target_path)"
    return target_path
end


"""
    loadWeightsAsDimArray(data::NCDataset, key_weights::String)

# Arguments:
- `data`: NCDataset containing weights, which have a single dimension
- `key_weights`: name of weights to load; 'wP' (performance weights), 'wI'
(independence weights), 'w' (overall weights)
"""
function loadWeightsAsDimArray(data::NCDataset, key_weights::String)
    src_name = dimnames(data[key_weights])[1]
    sources = Array(data[src_name])
    arr = DimArray(
        Array(data[key_weights]),
        (Dim{Symbol(src_name)}(sources)), metadata = Dict(data.attrib)
    )
    return arr
end

function loadWeightsFromJLD2(target_path::String)
    f = jldopen(target_path, "r")
    weights = f["weights"]
    close(f)
    return weights
end



""" 
    applyWeights(model_data::AbstractArray, weights::DimArray)

Compute the weighted average of model data 'data' with given weights 'weights'.
If weights were computed for a superset of the models in 'data', they are normalized
and applied to the subset. Only weights per model (not members) are considered
for now, in the future, members should be considered too.

# Arguments:
- `model_data`: model predictions. If given for model members, the predictions 
of each model are considered the average value of all members of the respective model.
- `weights`: if given for each member of a model, these will be summed up to 
yield one value per model.
"""
function applyWeights(model_data::AbstractArray, weights::DimArray)
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
        model_data = model_data[model = Where(x -> x in shared_models)]
        weights = weights[model = Where(x -> x in shared_models)]
        weights = weights ./ sum(weights)
        
        # sanity checks:
        if any(map(x -> !(x in shared_models), model_data.model))
            data_out = model_data[model = Where(x -> !(x in shared_models))]            
            models_out = dims(data_out, :model)
            @warn "No weights had been computed for $models_out"
        end 
        if length(shared_models) < length(models_weights)
            @warn "Weights were computed for a subset of the models of the given data. Weights were renormalized."
        end
    end
    return computeWeightedAvg(model_data; weights)
end


"""
    getModelLikelihoods(modelData::AbstractArray, distr::Distribution, diagnostic::String)

# Arguments:
- `modelData`:
- `distr`:
- `diagnostic`:
"""
function getModelLikelihoods(modelData::YAXArray, distr::Distribution, diagnostic::String)
    dim_symbol = hasdim(modelData, :model) ? :model : :member
    names_models = dims(modelData, dim_symbol)
    likelihoods = dim_symbol == :model ? 
        map(m -> Distributions.pdf(distr, DimArray(modelData)[model = At(m)]), names_models) :
        map(m -> Distributions.pdf(distr, DimArray(modelData)[member = At(m)]), names_models)

    variables = unique(filter(x->!ismissing(x), modelData.properties["variable_id"]))
    if length(variables) > 1
        throw(ArgumentError("Variable should be unique when computing likelihoods, given variables: $variables"))
    end
    return YAXArray(
        (Dim{dim_symbol}(Array(names_models)), Dim{:variable}([variables[1]]), Dim{:diagnostic}([diagnostic])),
        reshape(likelihoods, length(likelihoods), 1, 1),
        deepcopy(modelData.properties)
    )
end