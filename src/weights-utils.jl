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
    model_dim = Data.modelDim(data)
    equal_weighting = isnothing(weights)
    if equal_weighting
        weights = equalWeights(data; use_members = use_members_equal_weights)
    elseif !hasdim(weights, model_dim)
        msg = "weight vector must have same dimension as data (found: data: $model_dim, weights::$(dims(weights)))!"
        throw(ArgumentError(msg))
    end
    
    models_align = sort(collect(dims(data, model_dim))) == sort(collect(dims(weights, model_dim)))
    if !equal_weighting && !models_align
        weights, data = Data.alignWeightsAndData(data, weights)
    end
    weighted_data = @d data .* weights
    #weighted_avg = mapslices(x -> sum(skipmissing(x)), weighted_data; dims=(model_dim,))
    weighted_avg = mapslices(x -> sum(x), weighted_data; dims=(model_dim,))
    return weighted_avg
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
        if length(shared_models) < length(models)
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


# TODO: add config for independence weight, maybe abstract config struct, that can 
# be more general than ConfigWeights.
# return ClimWIP object
function climwipWeights(
    wP::MMEWeights, wI::MMEWeights; suffix_combined::String = ""
)
    wIP = wP.w .* wI.w
    wIP = wIP ./ sum(wIP)    
    # if isempty(name_combined)
    #     name_combined = join([wP.name, wI.name], "-")
    # end
    # if out == :only_combined
    #     return YAXArray(
    #         (Dim{:model}(lookup(wIP, :model)), Dim{:weight}([name_combined])), 
    #         reshape(wIP.data, :, 1)
    #     )
    # else
        # weights_arr = cat([wP.w, wI.w, wIP]..., dims=Dim{:weight}([wP.name, wI.name, name_combined]))
        # return weights_arr
    # end
    name_combined = isempty(suffix_combined) ? "wIP" : "wIP_" * suffix_combined
    weights_arr = cat([wP.w, wI.w, wIP]..., dims=Dim{:weight}([wP.name, wI.name, name_combined]))
    return ClimWIP(
        performance_distances = Dict(),
        independence_distances = Dict(),
        Di = YAXArray(rand(1,1)),
        Sij = YAXArray(rand(1,1)),
        w = weights_arr,
        config = ConfigWeights()
    )
end


function climwipWeights(
    weights_perform::AbstractVector{MMEWeights},
    wI::MMEWeights,
    suffix_combined_wP::String;
    suffix_combined::String = ""
)
    models = map(mme -> sort(lookup(mme.w, :model)), weights_perform) 
    if length(unique(models)) != 1
        throw(ArgumentError("Different performance weight vectors must be defined for the same models!"))
    end

    wP = similar(weights_perform[1].w)
    wP .= 1
    for mme_weights in weights_perform
        wP = @d wP .* mme_weights.w
    end
    wP = wP ./ sum(wP)
    name_wP = isempty(suffix_combined_wP) ? "wP" : "wP_" * suffix_combined_wP
    combined_wP = MMEWeights(w = wP, name = name_wP)
    return climwipWeights(combined_wP, wI; suffix_combined)
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


function performanceWeights(
    model_data::DataMap,
    obs_data::DataMap,
    diagnostics::AbstractVector{<:String};
    suffix_name::String = "",
    sigma_D::Number = 1/sqrt(2)
)
    dists_perform_all = Data.distancesData(model_data, obs_data, diagnostics)
    dists_perform_all = Data.summarizeMembers.(dists_perform_all)
    models = collect(lookup(dists_perform_all[1], :model))
    wPs = Vector{MMEWeights}(undef, length(diagnostics))
    for (i, d) in enumerate(diagnostics)
        wP_unnormalized = exp.(-Array(dists_perform_all[i].data)./(2*sigma_D^2))
        wP = wP_unnormalized ./ sum(wP_unnormalized)
        wPs[i] = MMEWeights(w = YAXArray((Dim{:model}(models),), wP), name=d * suffix_name)
    end
    return wPs
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


function brier(p_hat::Number, outcome::Int)
    if outcome != 1 && outcome != 0
        throw(ArgumentError("Event can either be observed (1) or not (1), input must be a number!"))
    end
    return (p_hat - outcome)^2
end


"""
    rps(p_hat::Vector{<:Number}, cat::Int)
"""
function rps(p_hat::Vector{<:Number}, cat::Int)
    if cat > length(p_hat)
        throw(ArgumentError("The length of the forecast probability vector ($p_hat) is smaller then the observed category ($cat)"))
    end
    cdf_p = cumsum(p_hat)
    obs = zeros(length(p_hat))
    obs[cat] = 1
    cdf_obs = cumsum(obs)
    return sum((cdf_p .- cdf_obs).^2)
end

"""
    crps(samples::AbstractVector{<:Union{Missing, Number}}, obs::Number)

Return Continuous Ranked Probability Score (CRPS) for predicted `samples` and observation `obs`.

Quantifies how good the forecast samples are a prediction for the given observation.
Use closed Form definition of CRPS.

# Arguments:
- `samples::AbstractVector{<:Union{Missing, Number}}`: vector of predictions.
- `obs::Number`:: single observation.
"""
function crps(samples::AbstractVector{<:Union{Missing, Number}}, obs::Number)
    # expected absolute error between forecast points and observation
    ev_obs = mean(abs.(samples .- obs))
    # expected absolute difference between two random forecast points (each sample pair)
    ev_pairs = mean(abs.(samples .- samples'))
    return ev_obs - 0.5 * ev_pairs
end


"""
    crps(samples::YAXArray{<:Union{Missing, Number}}, obs::YAXArray)

Compute Continuous Ranked Probability Score (CRPS) for predicted `samples` and observation `obs`.

# Arguments:
- `samples`: must have at least one more dimension than :model or :member (e.g. model x lon x lat)
- `obs`: must have the same dimensions as `samples` except for the model dimension (:model or :member)
which is only present in `samples`.
"""
function crps(samples::YAXArray{<:Union{Missing, Number}}, obs::YAXArray)
    model_dim = Data.modelDim(samples)
    dimensions = otherdims(samples, model_dim)
    if length(dimensions) == 0
        throw(ArgumentError("for observations that are single values, use single value as input arg to crps for obs! "))
    end
    if dimensions != dims(obs)
        throw(ArgumentError("Models and observations must be defined on same grid!"))
    end
    values = YAXArray(dimensions, ones(size(dimensions)))
    for (I, slice) in zip(CartesianIndices(obs), eachslice(samples; dims=dimensions))
        values[I] = crps(collect(slice), only(obs[I]))
    end
    return values
end


function crpss(
    samples_forecast::AbstractVector{<:Union{Missing, Number}}, 
    samples_baseline::AbstractVector{<:Union{Missing, Number}},
    obs::Number
)
    crps_fc = crps(samples_forecast, obs)
    crps_baseline = crps(samples_baseline, obs)
    return 1 - (crps_fc / crps_baseline)
end

function crpss(crps_fc::T, crps_baseline::T) where {T <: Union{Matrix, Number}}
    result = 1 .- (crps_fc ./ crps_baseline)
    return result
end

"""
   crpss(samples_forecast::YAXArray, samples_baseline::YAXArray, observations::YAXArray})
   
Compute Continuous Ranked Probability Skill Score for forecast with respect to baseline. 
    
For every timestep there are predicted samples of the forecast (`samples_forecast`) and of
the baseline `samples_baseline` (respectively one per model in the ensemble) and for every 
timestep there is exactly one observation. The CRPSS is based on the average CRPS across timesteps.

# Arguments:
- `samples_forecast::YAXArray`: must have dimension :time and :model or :member
- `samples_baseline::YAXArray`: must have dimension :time and :model or :member
- `observations::YAXArray`: must have only dimension :time
"""
function crpss(samples_forecast::YAXArray, samples_baseline::YAXArray, observations::YAXArray)
    Data.throwErrorIfDimMissing(samples_forecast, :time)
    Data.throwErrorIfDimMissing(samples_baseline, :time)
    Data.throwErrorIfDimMissing(observations, :time)
    timesteps = lookup(samples_forecast, :time)
    dims_fc = otherdims(samples_forecast, (:time, level_predictions))
    dims_bl = otherdims(samples_baseline, (:time, level_predictions))
    dims_obs = otherdims(observations, :time)
    if (dims_fc != dims_bl) || (dims_fc != dims_obs)
        throw(ArgumentError("forecast, baseline and observation must all have the same dimensions!"))
    end
    start_y, end_y = Dates.year.(timesteps[[1, end]])
    n_years = end_y - start_y + 1
    
    n_dims = length(dims_fc)
    crps_forecast = n_dims == 0 ? zeros(n_years) : zeros(size(dims_fc)..., n_years)
    crps_baseline = n_dims == 0 ? zeros(n_years) : zeros(size(dims_fc)..., n_years)
    indices = [Colon() for _ in 1:n_dims]
    for ts in range(1, n_years)
        if n_dims == 0
            crps_forecast[ts] = crps(samples_forecast[time = ts], observations[time=ts])
            crps_baseline[ts] = crps(samples_baseline[time = ts], observations[time=ts])      
        else
            crps_forecast[indices..., ts] = crps(samples_forecast[time = ts], observations[time=ts])
            crps_baseline[indices..., ts] = crps(samples_baseline[time = ts], observations[time=ts])
        end
    end
    # take mean over last dimension time
    crps_fc = dropdims(mean(crps_forecast, dims = n_dims + 1), dims = n_dims + 1) 
    crps_baseline = dropdims(mean(crps_baseline, dims = n_dims + 1), dims = n_dims + 1)
    result = crpss(crps_fc, crps_baseline)
    return n_dims > 0 ? YAXArray(dims_fc, result) : result
end


# Note: this function does not call crpss with the weighted samples because the weighted samples model dimension differs at every time step...!
"""
    crpss_weighted(samples::YAXArray, observations::YAXArray, weights::YAXArray; n::Int=1000)

Compute Continuous Ranked Probability Skill Score for weighted vs. unweighted predictions for timeseries.

For every timestep there are predicted samples, one per model in the ensemble and for 
every timestep there is exactly one observation. The CRPSS is based on the average CRPS
across timesteps.

# Arguments:
- `samples::YAXArray`: must have dimension :time and :model or :member
- `observations::YAXArray`: must have only dimension :time
- `weights::YAXArray`: must have same model dimension as `samples` (:model or :member)
- `n::Int=1000`: length sample vector for computing weighted samples
"""
function crpss_weighted(samples::YAXArray, observations::YAXArray, weights::YAXArray; n::Int=1000)
    level_predictions = Data.modelDim(samples)
    level_weights = Data.modelDim(weights)
    if level_weights != level_predictions
        throw(ArgumentError("Predictions and weights must be defined on the same level. Found: Predictions: $level_predictions, Weights: $level_weights."))
    end    
    Data.throwErrorIfDimMissing(samples, :time)
    Data.throwErrorIfDimMissing(observations, :time)
    timesteps = lookup(samples, :time)
    if Dates.year.(timesteps) != Dates.year.(lookup(observations, :time))
        throw(ArgumentError("forecast and observations must all be defined for the same time points!"))
    end
    start_y, end_y = Dates.year.(timesteps[[1, end]])
    n_years = end_y - start_y + 1
    dimensions = otherdims(samples, (:time, level_predictions))
    n_dims = length(dimensions)
    crps_weighted = n_dims == 0 ? zeros(n_years) : zeros(size(dimensions)..., n_years)
    crps_unweighted = n_dims == 0 ? zeros(n_years) : zeros(size(dimensions)..., n_years)
    indices = [Colon() for _ in 1:n_dims]
    for ts in range(1, n_years)
        forecast_samples = samples[time = ts]
        weighted_samples = weightSamples(forecast_samples, weights; n)
        #weighted_samples_indices = rand(Distributions.Categorical(weights), n)
        if n_dims == 0
            crps_weighted[ts] = crps(collect(weighted_samples), only(observations[time=ts]))
            crps_unweighted[ts] = crps(collect(forecast_samples), only(observations[time=ts]))
        else
            crps_weighted[indices..., ts] = crps(weighted_samples, observations[time=ts])
            crps_unweighted[indices..., ts] = crps(forecast_samples, observations[time=ts])
        end
    end
    # take mean over last dimension time
    crps_fc = dropdims(mean(crps_weighted, dims = n_dims + 1), dims = n_dims + 1) 
    crps_baseline = dropdims(mean(crps_unweighted, dims = n_dims + 1), dims = n_dims + 1)
    if n_dims > 0
        result = YAXArray(dimensions, crpss(crps_fc, crps_baseline))
    else
        result = crpss(only(crps_fc), only(crps_baseline))
    end
    return result
end

"""
    crpss_weighted(samples::YAXArray, observed::Number, weights::YAXArray)

Compute Continuous Ranked Probability Skill Score for weighted vs. unweighted predictions.

`samples` is one prediction per model in the ensemble, `observed` is the observed value.

# Arguments:
- `samples::YAXArray`: must have dimension :model or :member
- `observed::Number`:
- `weights::YAXArray`: must have same model dimension as `samples` (:model or :member)
- `n::Int=1000`: length sample vector for computing weighted samples
"""
function crpss_weighted(samples::YAXArray, observed::Number, weights::YAXArray; n::Int=1000)
    level_predictions = Data.modelDim(samples)
    level_weights = Data.modelDim(weights)
    if level_weights != level_predictions
        throw(ArgumentError("Predictions and weights must be defined on the same level. Found: Predictions: $level_predictions, Weights: $level_weights."))
    end
    weighted_samples = weightSamples(samples, weights; n)
    crps_weighted = crps(weighted_samples, observed)
    crps_unweighted = crps(samples, observed)
    return 1 - (crps_weighted / crps_unweighted)
end

"""
    weightSamples(samples::YAXArray, weights::YAXArray; n::Int = 1000)

Repeat each entry in `samples` proportional to respective value given in `weights` so that 
the returned vector has `n` entries.

# Arguments:
- `samples::YAXArray`: must have dimension :model or :member.
- `weights::YAXArray`: must have same model dimension as `samples` (:model or :member).
- `n::Int = 1000`: length of returned vector. 
"""
function weightSamples(samples::YAXArray, weights::YAXArray; n::Int = 1000)    
    level_predictions = Data.modelDim(samples)
    level_weights = Data.modelDim(weights)
    if level_weights != level_predictions
        throw(ArgumentError("Predictions and weights must be defined for the same models. Found: Predictions: $level_predictions, Weights: $level_weights."))
    end
    weights_n = rand(Distributions.Multinomial(n, weights))
    model_names = collect(lookup(samples, level_predictions))
    models = Vector{String}(undef, n)
    samples_mat = Array(samples)
    dimensions = otherdims(samples, level_predictions)
    n_dims = length(dimensions)
    weighted_samples = n_dims == 0 ? Array{eltype(samples_mat)}(undef, n) : 
        Array{eltype(samples_mat)}(undef, size(dimensions)..., n)
    indices = [Colon() for _ in 1:n_dims]
    pos = 1
    for (i, nb) in enumerate(weights_n)
        if n_dims == 0
            weighted_samples[pos : pos + nb - 1] .= samples_mat[indices..., i]
        else
            weighted_samples[indices..., pos : pos + nb - 1] .= samples_mat[indices..., i]
        end
        models[pos : pos + nb - 1] .= model_names[i]
        pos += nb
    end
    return YAXArray((dimensions..., Dim{level_predictions}(models)), weighted_samples)
end


"""
Compute weights proportional to area weighted mean squared error between model data and observations.
# Arguments:
- `data::YAXArray`: must have dimensions 'lon','lat', 'model'/'member' and 'diagnostic'.
- `obs::YAXArray`: must have dimensions 'lon','lat' and possibly 'model', and 'diagnostic'.
- `w_diagnostics::YAXArray`: must have dimension 'diagnostic with same values as in `model` and `data`;
values are normalized to sum to 1.
"""
function weightsAIC(data::YAXArray, obs::YAXArray, w_diagnostics::YAXArray)
    map(x -> Data.throwErrorIfDimMissing(x, [:diagnostic]), [data, obs, w_diagnostics])
    diagnostics = lookup(w_diagnostics, :diagnostic)
    n_a = (w_diagnostics./sum(w_diagnostics)) .* length(diagnostics)
    model_dim = Data.modelDim(data)
    n_models = length(lookup(data, model_dim))

    lls = zeros(n_models)
    for d in diagnostics
        mses = Data.distancesData(data[diagnostic = At(d)], obs[diagnostic = At(d)]; metric=:mse)
        mse_min = minimum(mses)
        lls .+= n_a[diagnostic = At(d)].data[] * log.(mse_min ./ mses)
    end
    likelihoods = exp.(lls)
    weights = likelihoods ./ sum(likelihoods)
    return YAXArray((dims(data, model_dim),), weights)
end


"""
Compute weights proportional to area weighted mean squared error between model data and observations.
# Arguments:
- `data::YAXArray`: must have dimensions 'lon','lat', 'model'/'member'.
- `obs::YAXArray`: must have dimensions 'lon','lat' and possibly 'model'.
"""
function weightsAIC(data::YAXArray, obs::YAXArray)
    mses = Data.distancesData(data, obs; metric=:mse)
    likelihoods = 1 ./ mses
    weights = likelihoods ./ sum(likelihoods)
    return YAXArray((dims(data, Data.modelDim(data)),), weights)
end


