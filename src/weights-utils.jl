"""
    weightedAvg(
        data::YAXArray; 
        weights::Union{YAXArray, Nothing} = nothing, 
        use_members_equal_weights::Bool = true
    )

Compute the average values for each (lon,lat) grid point in `data`, weighted
by `weights`. 

If no weight vector is provided, unweighted average is computed.


# Arguments:
- `data::YAXArray`: must have dimension 'model' or 'member'
- `weights::Union{YAXArray, Nothing} = nothing`: weights for models or individual members
- `use_members_equal_weights::Bool`:  if `weights` is nothing, all models receive 
equal weight. Then, if `use_members_equal_weights` is true, the number of members 
per model is taken into account, e.g. ["m1#run1", "m1#run2", "m1#run3", "m2#run1"] yields 
[1/4, 1/4, 1/4, 1/4] if `use_members_equal_weights = false` and if 
`use_members_equal_weights = true`, [1/2 * 1/3, 1/2 * 1/3, 1/2 * 3, 1/2] = [1/6, 1/6, 1/6, 1/2].
"""
function weightedAvg(
    data::YAXArray;
    weights::Union{YAXArray, Nothing} = nothing,
    use_members_equal_weights::Bool = true
)
    data = deepcopy(data)
    weights = deepcopy(weights)
    dim_symbol = Data.modelDim(data)
    equal_weighting = isnothing(weights)
    if equal_weighting
        weights = equalWeights(data; use_members = use_members_equal_weights)
    elseif !hasdim(weights, dim_symbol)
        msg = "weight vector must have same dimension as data (found: data: $dim_symbol, weights::$(dims(weights)))!"
        throw(ArgumentError(msg))
    end
    
    models_align = sort(collect(dims(data, dim_symbol))) == sort(collect(dims(weights, dim_symbol)))
    if !equal_weighting && !models_align
        weights, data = Data.alignWeightsAndData(data, weights)
    end
    weighted_data = @d data .* weights
    weighted_avg = sum(weighted_data, dims = dim_symbol)
    return Data.getAtModel(weighted_avg, dim_symbol, "combined")
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
        counts = Data.countMap(models)
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
    w_members = map(models) do m
        n_members = count(x -> startswith(x, m * MODEL_MEMBER_DELIM), members)
        w = only(weights[model=At(m)])
        fill(w / n_members, n_members)
    end
    return YAXArray((Dim{:member}(members),), reduce(vcat, w_members))
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
        weights = Data.summarizeMembersVector(weights; fn = sum)
    end
    if hasdim(model_data, :member)
        # take average over model predictions of members of same model
        model_data = Data.summarizeMembersVector(model_data; fn = Statistics.mean)
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
            @warn "Weights renormalized since computed for more models than models in given data."
        end
    end
    if hasdim(weights, :weight)
        results = map(x -> weightedAvg(model_data; weights=weights[weight=At(x)]), weights.weight)
        results_yax = cat(results...; dims=dims(weights, :weight))
    else
        results_yax = weightedAvg(model_data; weights)
    end
    return YAXArray(dims(results_yax), Array(results_yax), results_yax.properties)
end


function likelihoodWeights(
    vec::Vector{<:Number}, models::Vector{String}, distr::Distribution, name::String
)
    return likelihoodWeights(Distributions.pdf(distr, vec), models, name)
end

function likelihoodWeights(
    vec::Vector{<:Number}, models::Vector{String}, lh_fn::Function, name::String, args...; kwargs...
)
    likelihoods = map(x -> lh_fn(x, args...), vec)
    return likelihoodWeights(likelihoods, models, name)
end

function likelihoodWeights(
    likelihoods::Vector{<:Number}, models::Vector{String}, name::String
)
    weights = likelihoods ./ sum(likelihoods)
    return YAXArray((Dim{:model}(models), Dim{:weight}([name])), reshape(weights, :, 1))
end



function climwipWeights(
    wP::MMEWeights, wI::MMEWeights; out::Symbol = :only_combined, name_combined::String = ""
)
    w_combined = wP.w .* wI.w
    w_combined = w_combined ./ sum(w_combined)    
    if isempty(name_combined)
        name_combined = join([wP.name, wI.name], "-")
    end
    if out == :only_combined
        return YAXArray(
            (Dim{:model}(lookup(w_combined, :model)), Dim{:weight}([name_combined])), 
            reshape(w_combined.data, :, 1)
        )
    else
        weights_arr = cat([wP.w, wI.w, w_combined]..., dims=Dim{:weight}([wP.name, wI.name, name_combined]))
        return weights_arr
    end
end


function climwipWeights(
    weights_perform::AbstractVector{MMEWeights},
    wI::MMEWeights;
    out::Symbol = :only_combined,
    name_combined::String = ""
)
    models = map(mme -> sort(lookup(mme.w, :model)), weights_perform) 
    if length(unique(models)) != 1
        throw(ArgumentError("Different performance weight vectors must be defined for the same models!"))
    end

    wP = similar(weights_perform[1].w)
    wP .= 1
    name = ""
    for mme_weights in weights_perform
        wP = @d wP .* mme_weights.w
        name = join([name, mme_weights.name], "-")
    end
    wP = wP ./ sum(wP)
    combined_wP = MMEWeights(w=wP, name=name[2:end])
    return climwipWeights(combined_wP, wI; out, name_combined)
end

# functions for climwipWeights with historical performance:
"""
    climwipWeights(
        model_data::DataMap,
        obs_data::DataMap,
        config::ConfigWeights;
        suffix_name::String = "historical"
    )

Compute ClimwipWeights based on RMSE distances between observational data and model-model pairs.
"""
function climwipWeights(
    model_data::DataMap,
    obs_data::DataMap,
    config::ConfigWeights;
    suffix_name::String = "historical",
    norm_avg_members::Bool = true
)
    diagnostics_indep = Data.activeDiagnostics(config.independence)
    diagnostics_perform = Data.activeDiagnostics(config.performance)
    
    config = normalizeConfigWeight(config)
    weights_perform = Data.dict2YAX(config.performance)
    weights_indep = Data.dict2YAX(config.independence)

    dists_perform_all = Data.distancesData(model_data, obs_data, diagnostics_perform)
    dists_indep_all = Data.distancesModels(model_data, diagnostics_indep)

    Di = Data.generalizedDistances(dists_perform_all, diagnostics_perform, weights_perform; norm_avg_members)
    Sij = Data.generalizedDistances(dists_indep_all, diagnostics_indep, weights_indep; norm_avg_members)

    wP = performanceWeights(Di, config.sigma_performance)
    wI = independenceWeights(Sij, config.sigma_independence)

    w = wP .* wI
    w = w ./ sum(w)

    names_weights = String.(map(x -> join([x, suffix_name], "-"), ["wIP", "wI", "wP"]))
    weights_arr = cat([w, wI, wP]..., dims = Dim{:weight}(names_weights))
    
    # save dists and diagnostics as dictionaries
    dists_perform, dists_indep = Dict(), Dict()
    for (diag, dists) in zip(diagnostics_perform, dists_perform_all)
        dists_perform[diag] = dists
    end
    for (diag, dists) in zip(diagnostics_indep, dists_indep_all)
        dists_indep[diag] = dists
    end
    
    return ClimWIP(
        performance_distances = dists_perform,
        independence_distances = dists_indep,
        Di = Di,
        Sij = Sij,
        w = weights_arr,
        config = config
    )
end

"""
    climwipWeights(
        dists_indep_all::YAXArray, 
        dists_perform_all::YAXArray,
        config::ConfigWeights;
        suffix_name::Stirng = "climwip"
    )

Compute weight for each model in multi-model ensemble according to approach
from Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner,
Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming
from CMIP6 Projections When Weighting Models by Performance and
Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020):
995–1012. https://doi.org/10.5194/esd-11-995-2020. 

# Arguments:
- `dists_indep_all::YAXArray`: RMSEs between pairs of models for all 
combinations of variables and diagnostics; has dimensions 'member1', 'member2', 'variable', 
'diagnostic'.
- `dists_perform_all::YAXArray`: RMSEs between model and observational data for all 
combinations of variables and diagnostics; has dimensions 'member', 'variable', 'diagnostic'.
- `config::ConfigWeights`: Parameters specifiying the relative contributions of each 
combination of variable and diagnostic.
- `suffix_name::String`: added to name of each type of weights, s.t. names are: wP-suffix, 
wI-suffix, combined-suffix (default 'climwip').
"""
function climwipWeights(
    dists_indep_all::YAXArray,
    dists_perform_all::YAXArray,
    config::ConfigWeights;
    suffix_name::String = "historical"
)
    config = normalizeConfigWeight(config)
    weights_perform = Data.dict2YAX(config.performance)
    weights_indep = Data.dict2YAX(config.independence)

    Di = Data.generalizedDistances(dists_perform_all, weights_perform)
    Sij = Data.generalizedDistances(dists_indep_all, weights_indep)

    wP = performanceWeights(Di, config.sigma_performance)
    wI = independenceWeights(Sij, config.sigma_independence)

    w = wP .* wI
    w = w ./ sum(w)

    names = String.(map(x -> join([x, suffix_name], "-"), ["wIP", "wI", "wP"]))
    weights_arr = cat([w, wI, wP]..., dims = Dim{:weight}(names))
    return ClimWIP(
        performance_distances = dists_perform_all,
        independence_distances = dists_indep_all,
        Di = Di,
        Sij = Sij,
        w = weights_arr,
        config = config
    )
end


"""
    performanceWeights(Di::AbstractArray, sigmaD::Number)

Compute the ClimWIP performance weight for each model i with generalized distances `Di`.

```math
w^{P}_i = a_0 \\cdot \\frac{e^{-(\\frac{D_i}{\\sigma_D})^2}}{N}
```

# Arguments:
- `Di::YAXArray`: generalized distances Di for each model
- `sigmaD::Number`: free shape parameter for impact of performance weights
"""
function performanceWeights(Di::YAXArray, sigmaD::Number)
    performances = exp.(-(Di ./ sigmaD) .^ 2)
    return performances ./ sum(performances)
end

function performanceWeights(dists_perform_all::YAXArray, config::ConfigWeights)
    weights = Data.dict2YAX(Data.normalizeDict(config.performance))
    Di = Data.generalizedDistances(dists_perform_all, weights)
    return performanceWeights(Di, config.sigma_performance)
end


"""
    independenceWeights(Sij::YAXArray, sigmaS::Number)

Compute the independence weight for each model i with generalized distances between pairs of
models, `Sij`.

```math
w^{I}_i = a_0 \\cdot \\frac{1}{1 + \\sum_{j \\ne i} e^{-\\left( \\frac{S_{ij}}{\\sigma_S} \\right)^2}}
```

# Arguments:
- `Sij`: generalized distances for each pair of models with dimensions 'model1', 'model2'
- `sigmaS`: free shape parameter for impact of independence weights
"""
function independenceWeights(Sij::YAXArray, sigmaS::Number)
    # note: (+1 in eq. in paper is from when model is compared to itself since exp(0)=1)
    independences = sum(exp.(-(Sij ./ sigmaS) .^ 2), dims = :model2)
    independences = dropdims(independences, dims = :model2)
    independences = DimensionalData.set(independences, :model1 => :model)
    independences = 1 ./ independences
    return independences ./ sum(independences)
end

function independenceWeights(dists_indep_all::YAXArray, config::ConfigWeights)
    weights = Data.dict2YAX(Data.normalizeDict(config.independence))
    Sij = Data.generalizedDistances(dists_indep_all, weights)
    return independenceWeights(Sij, config.sigma_independence)
end


function normalizeConfigWeight(config::ConfigWeights)
    perform = Data.normalizeDict(config.performance)
    indep = Data.normalizeDict(config.independence)
    return ConfigWeights(perform, indep, config.sigma_performance, config.sigma_independence)
end