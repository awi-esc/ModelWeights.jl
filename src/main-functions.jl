import YAML
using DimensionalData
using Setfield


""" computeModelDataRMSE

# Arguments:
- `model_data`:
- `obs_data`: Observational data for computing performance weights.
- `config`:
- `ref_period_alias`:
"""
function computeModelDataRMSE(
    model_data::Vector{Data},
    obs_data::Vector{Data}, 
    config::ConfigWeights,
    ref_period_alias::String
)
    return computeDistancesAllDiagnostics(
        model_data, obs_data, config.performance, ref_period_alias, true
    )
end


"""
    computeWeights(
        model_data::Vector{Data}, dists_perform_all::DimArray, config::ConfigWeights
    )

Compute weight for each model in multi-model ensemble according to approach
from Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner,
Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming
from CMIP6 Projections When Weighting Models by Performance and
Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020):
995–1012. https://doi.org/10.5194/esd-11-995-2020. 

# Arguments:
- `model_data`: Models for data for computing independence and performance weights.
- `dists_perform_all`: distances for all variables and diagnostics, has dimensions, 
'member', 'variable', 'diagnostic'.
- `config`: parameters specifiying the relative contributions of each 
combination of variable and diagnostic.
"""
function computeWeights(
    model_data::Vector{Data}, dists_perform_all::DimArray, config::ConfigWeights
)
    weights_perform = normalizeWeightsVariables(config.performance)  
    weights_indep = normalizeWeightsVariables(config.independence)
    #keys_weights_perform = allcombinations(dims(weights_perform, :variable), dims(weights_perform, :diagnostic))
    keys_weights_indep = allcombinations(dims(weights_indep, :variable), dims(weights_indep, :diagnostic))
    ref_period_alias, ref_period_timerange = getRefPeriodAsTimerangeAndAlias(
        map(x -> x.meta.attrib, model_data), config.ref_perform_weights
    )
    # sanity checks for input arguments
    msg(x) =  "For computation of $x weights: Make sure that data is provided 
    for the given reference period ($(config.ref_perform_weights)) and combination of 
        variables and diagnostic for which (balance) weights were specified!"
    # if !isValidDataAndWeightInput(model_data, keys_weights_perform, ref_period_alias)
    #     throw(ArgumentError(msg("performance")))
    # end
    # if !isValidDataAndWeightInput(obs_data, keys_weights_perform, ref_period_alias)
    #     throw(ArgumentError(msg("performance")))
    # end
    if !isValidDataAndWeightInput(model_data, keys_weights_indep, ref_period_alias)
        throw(ArgumentError(msg("independence")))
    end

    dists_indep_all = computeDistancesAllDiagnostics(
        model_data, nothing, config.independence, config.ref_indep_weights, false
    )
    Di = computeGeneralizedDistances(dists_perform_all, weights_perform, true)
    Sij = computeGeneralizedDistances(dists_indep_all, weights_indep, false)

    performances = performanceParts(Di, config.sigma_performance)
    independences = independenceParts(Sij, config.sigma_independence)
    weights = performances ./ independences
    weights = weights ./ sum(weights)
    setRefPeriodInWeightsMetadata!(weights.metadata, ref_period_alias, ref_period_timerange)
    
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
    setRefPeriodInWeightsMetadata!(wP.metadata, ref_period_alias, ref_period_timerange)
    
    if !isempty(config.target_path)
        target_path = validateConfigTargetPath(config.target_path)
        # update config accordingly since it is also saved together with the weights
        config = @set config.target_path = target_path
    end
    
    w_members = distributeWeightsAcrossMembers(weights);
    model_weights =  Weights(
        performance_distances = dists_perform_all,
        independence_distances = dists_indep_all, 
        Di = Di,
        Sij = Sij,
        wP = wP,
        wI = wI,
        w =  weights,
        w_members = w_members,
        config = config
        # overall = (wP .* wI) ./ sum(wP .* wI), # just for sanity check
    )
    @debug weights
    if !isempty(config.target_path)
        filename = basename(target_path)
        if endswith(filename, ".nc")
            saveWeightsAsNCFile(model_weights, target_path)
        elseif endswith(filename, ".jld2")
            saveWeightsAsJuliaObj(model_weights, target_path)
        else
            @warn "Weights not saved - they can only be saved as .nc or .jld2 files!"
        end 
    end
    return model_weights
end


"""
    loadDataFromESMValToolConfigs(
        path_data::String,
        path_recipes::String;
        dir_per_var::Bool=true,
        is_model_data::Bool=true,
        subset::Union{Dict, Nothing}=nothing
    )

Loads the data from the config files (ESMValTool recipes) located at 'path_recipes'.
For each variable, experiment, statistic and timerange (alias) a different
DimArray is loaded.

# Arguments:
- `path_data`:  path to location where preprocessed data is stored; if 
dir_per_var is true, paths to directories that contain one or
more subdirectories that each contains a directory 'preproc' with the
preprocessed data. If dir_per_var is false, path_data is path to directory that
 directly contains the 'preproc' subdirectory.
- `path_recipes`: path to directory that contains one or ESMValTool recipes 
used as config files
- `dir_per_var`: if true, directory at path_data has subdirectories, one for
each variable (they must end with _ and the name of the variable), otherwise
data_path points to the directory that has a subdirectory 'preproc'.
- `is_model_data`: set true for model data, false for observational data
- `subset`: specifies the properties of the subset of data to be loaded. 
The following keys are considered:  `models` (used to load only specific set of models 
or members of models), `projects` (used to load only data from a given set of
projects, e.g. for loading only CMIP5-data), `timeranges` and `aliases` 
(super set is loaded, i.e. data that corresponds to either a given timerange or
a given alias will be loaded), `variables`, `statistics`, `subdirs` and `level_shared_models` 
(if set to MEMBER/MODEL only data loaded from model members/models shared across all experiments and variables is loaded).
- `preview`: used to pre-check which data will be loaded before actually loading
it. 
"""
function loadDataFromESMValToolConfigs(
    path_data::String,
    path_recipes::String;
    dir_per_var::Bool=true,
    is_model_data::Bool=true,
    subset::Union{Dict, Nothing}=nothing,
    preview::Bool=false
)
    attributes = getMetaAttributesFromESMValToolConfigs(path_recipes; constraint=subset)
    meta_data = Dict{String, MetaData}()
    for attrib in attributes
        meta = buildMetaData(
            attrib, path_data, dir_per_var, is_model_data; constraint=subset
        )
        addMetaData!(meta_data, meta)
    end
    if !isnothing(subset) && !isnothing(get(subset, "level_shared_models", nothing))
        reduceMetaDataSharedModels!(meta_data, subset["level_shared_models"])
    end
    output = preview ? meta_data : loadDataFromMetadata(meta_data, is_model_data)
    return output
end


"""
    loadDataFromYAML(
        path_config::String;
        is_model_data::Bool=true,
        subset::Union{Dict, Nothing}=nothing,
        preview::Bool=false
    )

Load a `DataMap`-instance with data as specified in yaml file at `path_config`.

# Arguments:
- `path_config`: path to yaml config file
- `is_model_data`: true for model data (default), false for observational data
- `subset`: specifies the properties of the subset of data to be loaded. These 
have to apply to all loaded datasets specified in the config yaml file. 
The following keys are considered:  `models` (used to load only specific set of models 
or members of models), `projects` (used to load only data from a given set of
projects, e.g. for loading only CMIP5-data), `timeranges` and `aliases` 
(super set is loaded, i.e. data that corresponds to either a given timerange or
a given alias will be loaded), `variables`, `statistics`, `subdirs`, `level_shared_models` 
(if set to MEMBER/MODEL only data loaded from model members/models shared across all
experiments and variables is loaded) and `dir_per_var`.
- `preview`: used to pre-check which data will be loaded before actually loading
it. 
"""
function loadDataFromYAML(
    path_config::String;
    is_model_data::Bool=true,
    subset::Union{Dict, Nothing}=nothing,
    preview::Bool=false
)
    meta_data = getMetaDataFromYAML(path_config, is_model_data; arg_constraint = subset)
    if !isnothing(subset) && !isnothing(get(subset, "level_shared_models", nothing))
        reduceMetaDataSharedModels!(meta_data, subset["level_shared_models"])
    end
    output = preview ? meta_data : loadDataFromMetadata(meta_data, is_model_data)
    return output
end


