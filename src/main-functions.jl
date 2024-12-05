import YAML
using DimensionalData


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
    model_data::Data, obs_data::Data, config::ConfigWeights
)
    weights_perform = normalizeWeightsVariables(config.performance)  
    weights_indep = normalizeWeightsVariables(config.independence)
    
    # sanity checks for input arguments
    keys_weights_perform = allcombinations(dims(weights_perform, :variable), dims(weights_perform, :diagnostic))
    keys_weights_indep = allcombinations(dims(weights_indep, :variable), dims(weights_indep, :diagnostic))
    
    ref_period_alias, ref_period_timerange = getRefPeriodAsTimerangeAndAlias(
        map(x -> x.attrib, model_data.meta), config.ref_period
    )
    
    msg =  x -> "For computation of $x weights: Make sure that data is provided 
    for the given reference period ($(config.ref_period)) and combination of 
        variables and diagnostic for which (balance) weights were specified!"
    if !isValidDataAndWeightInput(model_data, keys_weights_perform, ref_period_alias)
        throw(ArgumentError(msg("performance")))
    end
    if !isValidDataAndWeightInput(obs_data, keys_weights_perform, ref_period_alias)
        throw(ArgumentError(msg("performance")))
    end
    if !isValidDataAndWeightInput(model_data, keys_weights_indep, ref_period_alias)
        throw(ArgumentError(msg("independence")))
    end

    dists_perform_all = computeDistancesAllDiagnostics(
        model_data, obs_data, config.performance, ref_period_alias, true
    )
    dists_indep_all = computeDistancesAllDiagnostics(
        model_data, nothing, config.independence, ref_period_alias, false
    )
    Di = computeGeneralizedDistances(dists_perform_all, weights_perform, true)
    Sij = computeGeneralizedDistances(dists_indep_all, weights_indep, false)

    performances = performanceParts(Di, config.sigma_performance)
    independences = independenceParts(Sij, config.sigma_independence)
    weights = performances ./ independences;
    weights = weights ./ sum(weights);
    setRefPeriodInWeightsMetadata!(weights.metadata, ref_period_alias, ref_period_timerange)
    
    wI = independences ./ sum(independences)
    wP = performances ./ sum(performances)
    setRefPeriodInWeightsMetadata!(wP.metadata, ref_period_alias, ref_period_timerange)
    
    climwip_weights =  ClimwipWeights(
        performance_distances = dists_perform_all,
        independence_distances = dists_indep_all, 
        Di = Di,
        Sij = Sij,
        wP = wP,
        wI = wI,
        w =  weights
        #overall = (wP./wI)./sum(wP./wI), # just for sanity check
        )
    logWeights(weights.metadata);
    if !isempty(config.target_dir)
        saveWeights(climwip_weights, config.target_dir)
    end
    return climwip_weights
end


"""
    loadDataFromESMValToolConfigs(
        path_data::String,
        path_recipes::String;
        dir_per_var::Bool=true,
        is_model_data::Bool=true,
        only_shared_models::Bool=false,
        subset::Union{Dict{String, Vector{String}}, Nothing}=nothing
    )

Loads the data from the config files located at 'config_paths'. For necessary
structure of config files, see TODO. For each variable, experiment, statistic
and timerange (alias) a different DimArray is loaded.

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
- `only_shared_models`: if true only data loaded from model members shared 
across all experiments and variables.
- `subset`: dictionary specifying the subset of data to be loaded, ech key maps 
to a vector of Strings; the following keys are handled: 
'variable', 'statistic', 'alias', `timerange`, `exp` (fields of struct 
`MetaAttrib`), `models` which maps to Vector of specific models or members of
models (default is empty) and `projects` which maps to Vector of Strings, 
for instance ["CMIP6"] to load only CMIP6-data.
"""
function loadDataFromESMValToolConfigs(
    path_data::String,
    path_recipes::String;
    dir_per_var::Bool=true,
    is_model_data::Bool=true,
    only_shared_models::Bool=false,
    subset::Union{Dict{String, Vector{String}}, Nothing}=nothing
)
    attributes = getMetaAttributesFromESMValToolConfigs(path_recipes; subset)
    meta_data = buildMetaData(
        attributes, path_data, dir_per_var, is_model_data; subset
    )
    data = loadDataFromMetadata(meta_data, is_model_data, only_shared_models)
    return data
end


function loadDataFromYAML(
    path_config::String;
    dir_per_var::Bool=true,
    is_model_data::Bool=true,
    only_shared_models::Bool=false,
    subset::Union{Dict{String, Vector{String}}, Nothing}=nothing
)
    meta_data = getMetaDataFromYAML(path_config, dir_per_var, is_model_data; subset)
    data = loadDataFromMetadata(meta_data, is_model_data, only_shared_models)
    return data
end
