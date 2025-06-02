import YAML

using DimensionalData
using Serialization
using Setfield


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
    weights_perform = normalizeWeightsVariables(config.performance)
    weights_indep = normalizeWeightsVariables(config.independence)

    Di = computeGeneralizedDistances(dists_perform_all, weights_perform)
    Sij = computeGeneralizedDistances(dists_indep_all, weights_indep)

    performances = performanceParts(Di, config.sigma_performance)
    independences = independenceParts(Sij, config.sigma_independence)
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
    dists_perform = computeDistancesAllDiagnostics(data.models, data.obs, config_weights.performance)
    dists_indep = computeDistancesAllDiagnostics(data.models, config_weights.independence)
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


"""
    loadDataFromESMValToolRecipes(
        path_data::String,
        path_recipes::String;
        dir_per_var::Bool=true,
        is_model_data::Bool=true,
        subset::Union{Dict, Nothing}=nothing,
        preview::Bool=false
    )

Loads the data from the config files (ESMValTool recipes) located at 'path_recipes'.
For each variable, experiment, statistic and timerange (alias), an instance of `Data`
is created.

# Arguments:
- `path_data`:  path to location where preprocessed data is stored; if 
dir_per_var is true, paths to directories that contain one or
more subdirectories that each contains a directory 'preproc' with the
preprocessed data. If dir_per_var is false, path_data is path to directory that
 directly contains the 'preproc' subdirectory.
- `path_recipes`: path to directory that contains one or ESMValTool recipes 
used as config files.
- `dir_per_var`: if true (default), directory at path_data has subdirectories, one for
each variable (they must end with _ and the name of the variable), otherwise
data_path points to the directory that has a subdirectory 'preproc'.
- `is_model_data`: set true for model data, false for observational data
- `subset`: specifies the properties of the subset of data to be loaded. 
The following keys are considered:  `models` (used to load only specific set of models 
or members of models), `projects` (used to load only data from a given set of
projects, e.g. for loading only CMIP5-data), `timeranges` and `aliases`.
(super set is loaded, i.e. data that corresponds to either a given timerange or
a given alias will be loaded), `variables`, `statistics`, `subdirs` and `subset_shared` .
(if set to MEMBER/MODEL only data loaded from model members/models shared across all experiments and variables is loaded).
- `preview`: used to pre-check which data will be loaded before actually loading
it. 
"""
function loadDataFromESMValToolRecipes(
    path_data::String,
    path_recipes::String;
    dir_per_var::Bool = true,
    is_model_data::Bool = true,
    subset::Union{Dict,Nothing} = nothing,
    preview::Bool = false,
)
    attributes = getMetaAttributesFromESMValToolConfigs(path_recipes; constraint = subset)
    meta_data = Dict{String, Dict{String, Any}}()
    for meta in attributes
        meta = Dict{String, Any}(meta)
        meta["_paths"] = getPathsToData(
            meta,
            path_data,
            dir_per_var,
            is_model_data;
            constraint = subset,
        )
        if !isempty(meta["_paths"])
            id = meta["_id"]
            if haskey(meta_data, id)
                meta_data[id]["_paths"] = mergeMetaDataPaths(meta_data[id], meta)
            else
                meta_data[id] = meta
            end
        end
    end
    if !isnothing(subset) && !isnothing(get(subset, "subset_shared", nothing))
        filterPathsSharedModels!(meta_data, subset["subset_shared"])
    end
    return preview ? meta_data : loadDataFromMetadata(meta_data, is_model_data)
end


"""
    loadDataFromYAML(
        path_config::String;
        is_model_data::Bool = true,
        subset::Union{Dict,Nothing} = nothing,
        preview::Bool = false
    )

Load a `DataMap`-instance that contains the data specified in yaml file at `path_config`, 
potentially constraint by values in `subset`.

# Arguments:
- `preview::Bool`: if true (default: false), return metadata without actually 
loading any data.
"""
function loadDataFromYAML(
    path_config::String;
    is_model_data::Bool = true,
    subset::Union{Dict,Nothing} = nothing,
    preview::Bool = false
)
    return loadDataFromYAML(YAML.load_file(path_config); is_model_data, subset, preview)
end


function loadDataFromYAML(
    content::Dict;
    is_model_data::Bool = true,
    subset::Union{Dict,Nothing} = nothing,
    preview::Bool = false,
)
    meta_data = getMetaDataFromYAML(content, is_model_data; arg_constraint = subset)
    if isempty(meta_data)
        msg = "No metadata found for subset: $subset, $content (model data: $is_model_data)!"
        throw(ArgumentError(msg))
    end
    return preview ? meta_data : loadDataFromMetadata(meta_data, is_model_data)
end