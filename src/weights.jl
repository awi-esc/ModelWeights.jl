"""
    performance(generalizedDistances::AbstractArray, sigmaD::Number)

Compute the performance part of the overall weight for each model i defined as:

```math
w^{P}_i = \\frac{e^{-(\\frac{D_i}{\\sigma_D})^2}}{N}
```

# Arguments:
- `generalizedDistances::AbstractArray`: contains generalized distances Di for each model
- `sigmaD::Number`: free shape parameter for impact of performance weights
"""
function performance(generalizedDistances::AbstractArray, sigmaD::Number)
    return exp.(-(generalizedDistances ./ sigmaD) .^ 2)
end


"""
    independence(generalizedDistances::AbstractArray, sigmaS::Number)

Compute the independence part of the overall weight for each model i defined as:

```math
w^{I}_i = \\frac{1}{1 + \\sum_{j \\ne i} e^{-\\left( \\frac{S_{ij}}{\\sigma_S} \\right)^2}}
```

# Arguments:
- `generalizedDistances`: contains generalized distances S_{i,j} for each model 
(i.e. not model members) pair has two dimensions, 'model1' and 'model2'
- `sigmaS`: free shape parameter for impact of independence weights
"""
function independence(generalizedDistances::AbstractArray, sigmaS::Number)
    # note: (+1 in eq. in paper is from when model is compared to itself since exp(0)=1)
    indep_parts = sum(exp.(-(generalizedDistances ./ sigmaS) .^ 2), dims = :model2)
    indep_parts = dropdims(indep_parts, dims = :model2)
    indep_parts = DimensionalData.set(indep_parts, :model1 => :model)
    return indep_parts
end


# TODO: check use_members_equal_weights not only taken into account for equal weights?! Also for weights on level of models?
"""
    weightedAvg(
        data::YAXArray; 
        weights::Union{DimArray, Nothing}=nothing, 
        use_members_equal_weights::Bool=true
    )

Compute the average values for each (lon,lat) grid point in `data`, weighted
by `weights`. 

If no weight vector is provided, unweighted average is computed.


# Arguments:
- `data::YAXArray`: must have dimension 'model' or 'member'
- `weights::Union{YAXArray, Nothing}=nothing`: weights for models or individual members
- `use_members_equal_weights::Bool`:  if `weights` is nothing, all models receive 
equal weight. If `use_members_equal_weights` is true, the number of members 
per model is taken into account, s.t. each model receives equal weight, which
is distributed among the respective members.
"""
function weightedAvg(
    data::YAXArray;
    weights::Union{YAXArray, Nothing} = nothing,
    use_members_equal_weights::Bool = true
)
    da = DimArray(Array{Union{Missing,Float64}}(undef, size(data)), dims(data))
    da .= Array(data)
    dim_symbol = hasdim(data, :member) ? :member : :model
    models_data = collect(dims(data, dim_symbol))

    equal_weights = isnothing(weights)
    if equal_weights
        weights = equalWeights(data; use_members = use_members_equal_weights)
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
            msg = "No weights were computed for follwoing models, thus not considered in the weighted average:"
            @warn msg data_no_weights
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
    @assert isapprox(sum(weights), 1; atol = 10^-4)
    # there shouldnt be data for which there is no weight since it was filtered out above
    @assert sort(Array(models_weights)) == sort(Array(dims(data, dim_symbol)))

    # readjust weights for missing values, if one model has a missing value, 
    # it is ignored in the weighted average for that particular lon,lat-position.
    not_missing_vals = mapslices(x -> (ismissing.(x) .== 0), da; dims = (dim_symbol,))
    w_temp = mapslices(
        x -> (x .* weights) ./ sum(x .* weights),
        not_missing_vals,
        dims = (dim_symbol,),
    )
    w_temp = replace(w_temp, NaN => missing)

    for m in models_weights
        value = Data.getAtModel(da, dim_symbol, m) .* Data.getAtModel(w_temp, dim_symbol, m)
        Data.putAtModel!(da, dim_symbol, m, value)
    end
    weighted_avg = mapslices(x -> sum(skipmissing(x)), da, dims = (dim_symbol,))
    weighted_avg = dropdims(weighted_avg, dims = dim_symbol)

    # set to missing when value was missing for ALL models
    n_models = length(dims(da, dim_symbol))
    all_missing = dropdims(
        mapslices(x -> sum(ismissing.(x)) == n_models, da, dims = (dim_symbol,)),
        dims = dim_symbol,
    )
    if !any(all_missing)
        result = weighted_avg
    else
        result = Array{Union{Missing,Float64}}(undef, size(weighted_avg))
        s = repeat([:], length(size(weighted_avg)))
        result[s...] = Array(weighted_avg)
        result[all_missing] .= missing
    end
    meta = deepcopy(data.properties)
    meta["_statistic"] = equal_weights ? "unweighted-avg" : "weighted-avg"
    meta["_id"] = Data.buildMetaDataID(meta)
    return YAXArray(dims(weighted_avg), result, meta)
end


"""
    equalWeights(data::YAXArray; use_members::Bool=true)

Create a weight vector, with equal weight for each MODEL. Distribute weight across
model members if dimension=:member and use_members is true. If use_member is false, 
each model member is considered as standalone model and all receive the same weight.
"""
function equalWeights(data::YAXArray; use_members::Bool = true)
    dimension = hasdim(data, :member) ? :member : :model
    dimnames = dims(data, dimension)
    models = collect(map(x -> split(x, MODEL_MEMBER_DELIM)[1], dimnames))
    n_models = length(unique(models))
    if dimension == :member
        # make sure that the number of members per model is considered
        counts = countmap(models)
        n_members = length(dimnames)
        w =
            use_members ? [1 / n_models * 1 / counts[m] for m in models] :
            repeat([1 / n_members], n_members)
    elseif dimension == :model
        w = [1 / n_models for _ in range(1, n_models)]
    else
        throw(ArgumentError("To compute weights data must have dim :model or dim :member"))
    end
    return YAXArray((dims(data, dimension),), w)
end


"""
    distributeWeightsAcrossMembers(weights::YAXArray, members::Vector{String})

Equally distribute weight for each model over its members.
"""
function distributeWeightsAcrossMembers(weights::YAXArray, members::Vector{String})
    models = Array(dims(weights, :model))
    counts_models = map(
        m -> length(findall(x -> startswith(x, m * MODEL_MEMBER_DELIM), members)),
        models,
    )
    w_members = []
    for (i, m) in enumerate(models)
        n_members = counts_models[i]
        w = only(weights[model=At(m)])
        for _ in range(1, n_members)
            push!(w_members, w / n_members)
        end
    end
    return YAXArray((Dim{:member}(members),), w_members)
end


"""
    saveWeightsAsNCFile(weights::ClimWIP; target_path::String)
TODO: needs to be updated, deprecated
# Arguments:
- `weights`: ClimWIP object to be saved.
- `target_path`: Path to where weights shall be stored.
"""
function saveWeightsAsNCFile(weights::ClimWIP, target_path::String)
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
    for name in fieldnames(ClimWIP)
        if String(name) in ["w", "wP", "wI", "Di"]
            v = defVar(ds, String(name), Float64, ("model",))
            data = getfield(weights, name)
            v[:] = Array(data)
        elseif String(name) == "independence_distances"
            v = defVar(
                ds,
                String(name),
                Float64,
                ("member", "member", "variable", "diagnostic"),
            )
            v[:, :, :, :] = Array(weights.independence_distances)
        elseif String(name) == "Sij"
            v = defVar(ds, String(name), Float64, ("model", "model"))
            v[:, :] = Array(weights.Sij)
        elseif String(name) == "performance_distances"
            v = defVar(ds, String(name), Float64, ("member", "variable", "diagnostic"))
            v[:, :, :] = Array(weights.performance_distances)
        elseif String(name) == "config"
            # configuration is added as global attributes:
            for config_name in fieldnames(ConfigWeights)
                config_val = getfield(weights.config, config_name)
                config_key = String(config_name)
                if config_key in ["performance", "independence"]
                    for (k, v) in config_val
                        ds.attrib["w_"*config_key*"_"*k] = v
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
    return target_path
end


function writeWeightsToDisk(weights::ClimWIP, target_path::String)
    target_path = individuatePath(target_path)
    config = weights.config
    config = @set config.target_path = target_path
    weights = @set weights.config = config
    if endswith(target_path, ".jld2")
        jldsave(target_path; weights = weights)
    elseif endswith(target_path, ".nc")
        saveWeightsAsNCFile(model_weights, target_path)
    else
        serialize(target_path, weights)
    end
    @info "saved weights to: $(target_path)"
    return target_path
end


# TODO check if still necessary
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
        (Dim{Symbol(src_name)}(sources)),
        metadata = Dict(data.attrib),
    )
    return arr
end


""" 
    applyWeights(model_data::YAXArray, weights::YAXArray)

Compute the weighted average of model data `data` with given `weights`.

If weights were computed for a superset of the models in `data`, they are normalized
and applied to the subset. Only weights per model (not members) are considered
for now, in the future, members should be considered too.

# Arguments:
- `model_data::YAXArray`: model predictions. If given for model members, the predictions 
of each model are considered the average value of all members of the respective model.
- `weights::YAXArray`: if given for each member of a model, these will be summed up to 
yield one value per model.
"""
function applyWeights(model_data::YAXArray, weights::YAXArray)
    if hasdim(weights, :member)
        weights = Data.summarizeEnsembleMembersVector(weights, true; fn = sum)
    end
    if hasdim(model_data, :member)
        # take average over model predictions of members of same model
        model_data = Data.summarizeEnsembleMembersVector(model_data, true; fn = Statistics.mean)
    end
    models_weights = dims(weights, :model)
    models = dims(model_data, :model)
    shared_models = intersect(models_weights, models)

    if length(shared_models) == 0
        throw(ArgumentError("No models given for which weights had been computed!"))
    else
        model_data = model_data[model = Where(x -> x in shared_models)]
        weights = weights[model = Where(x -> x in shared_models)]
        weights = weights ./ sum(weights; dims=:model)
        if length(shared_models) < length(models_weights)
            @warn "Weights renormalized since computed for subset of models in given data."
        end
    end
    results = map(x -> weightedAvg(model_data; weights=weights[weight=At(x)]), weights.weight)
    results_yax = cat(results...; dims=dims(weights, :weight))
    return YAXArray(dims(results_yax), Array(results_yax), results_yax.properties)
end


function likelihoodWeights(
    vec::Vector{<:Number}, models::Vector{String}, distr::Distribution, name::String
)
    likelihoods = Distributions.pdf(distr, vec)
    weights = likelihoods ./ sum(likelihoods)
    return YAXArray((Dim{:model}(models), Dim{:weight}([name])), reshape(weights, :, 1))
end


"""
    climwipWeights(
        dists_indep_all::YAXArray, 
        dists_perform_all::YAXArray,
        config::ConfigWeights
    )

Compute weight for each model in multi-model ensemble according to approach
from Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner,
Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming
from CMIP6 Projections When Weighting Models by Performance and
Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020):
995–1012. https://doi.org/10.5194/esd-11-995-2020. 

# Arguments:
- `dists_indep_all::YAXArray`: RMSEs between pairs of models for all 
combinations of variables and diagnostics; has dimensions 'member1', 'member2', 
'variable', 'diagnostic'.
- `dists_perform_all::YAXArray`: RMSEs between model and observational 
data for all combinations variables and diagnostics; has dimensions 'member', 
'variable', 'diagnostic'.
- `config::ConfigWeights`: Parameters specifiying the relative contributions 
of each combination of variable and diagnostic.
- `suffix::String`: added to name of each type of weights, s.t. names are: 
wP-suffix, wI-suffix, combined-suffix.
"""
function climwipWeights(
    dists_indep_all::YAXArray,
    dists_perform_all::YAXArray,
    config::ConfigWeights,
    suffix::String
)
    weights_perform = Data.normalize(config.performance)
    weights_indep = Data.normalize(config.independence)

    Di = Data.generalizedDistances(dists_perform_all, weights_perform)
    Sij = Data.generalizedDistances(dists_indep_all, weights_indep)

    performances = performance(Di, config.sigma_performance)
    independences = independence(Sij, config.sigma_independence)
    weights = performances ./ independences
    weights = weights ./ sum(weights)
    # consider performance and independence weights independently
    # for performance weights, we assume that all models have the same degree of dependence
    # among each other (e.g. all are compeletely independent), i.e. we can 
    # just consider the performance Parts (the denominator would be the same for all models)
    norm_p = sum(performances)
    wP = performances ./ norm_p

    # for independence weights, we assume that all models perform equally well, i.e. 
    # Di = Dj for all models i, j. Thus, the numerator would be the same for all models, 
    # we just set Di=0 for all models, i.e. the numerator is 1 for all models
    norm_i = sum(1 ./ independences)
    wI = (1 ./ independences) ./ norm_i

    if hasdim(dists_perform_all, :member)
        members = Array(dims(dists_perform_all, :member))
        w_members = distributeWeightsAcrossMembers(weights, members)
    else
        w_members = weights
    end
    names = String.(map(x -> join([x, suffix], "-"), ["wP", "wI", "combined"]))
    weights_arr = cat([wP, wI, weights]..., dims=Dim{:weight}(names))
    model_weights = ClimWIP(
        performance_distances = dists_perform_all,
        independence_distances = dists_indep_all,
        Di = Di,
        Sij = Sij,
        w = weights_arr,
        w_members = w_members,
        config = config,
    )
    return model_weights
end


function climwipWeights(
    data::ClimateData, config_weights::ConfigWeights; suffix::String="climwip"
)
    dists_perform = Data.distancesData(data.models, data.obs, config_weights.performance)
    dists_indep = Data.distancesModels(data.models, config_weights.independence)
    return climwipWeights(dists_indep, dists_perform, config_weights, suffix)
end


function addClimwipWeights!(
    data::ClimateData, config_weights::ConfigWeights; suffix::String="climwip"
)
    weights = climwipWeights(data, config_weights; suffix)
    addWeights!(data, weights.w)
    return nothing
end


function addClimwipWeights!(data::ClimateData, weights::ClimWIP)
    addWeights!(data, weights.w)
    return nothing
end



"""
    addLikelihoodWeights!(
        data::ClimateData, 
        values::Vector{<:Number}, 
        models::Vector{String}, 
        distr::Distribution,
        id::String
    )

Compute the likelihood weights for `values` with respect to the distribution `distr` and 
add computed weights to weights Array of `data` with name `id`. 
If `id` is already present, nothing is added.

"""
function addLikelihoodWeights!(
    data::ClimateData, 
    values::Vector{<:Number}, 
    models::Vector{String}, 
    distr::Distribution,
    id::String
)
    addWeights!(data, likelihoodWeights(values, models, distr, id))
    return nothing
end



function addWeights!(data::ClimateData, weights::YAXArray)
    if isnothing(data.weights)
        data.weights = weights
    else
        names_weights = Array(data.weights.weight)
        names_new_weights = Array(weights.weight)
        isPresentID = map(x -> x in names_weights, names_new_weights)
        if any(isPresentID)
            @warn "$(names_new_weights[isPresentID]) already present in weights! Nothing added."
            return nothing
        end
        lookups = vcat(names_weights, Array(weights.weight))
        data.weights = cat(data.weights, weights; dims=Dim{:weight}(lookups))
    end
    return nothing
end

