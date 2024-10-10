import YAML
using DimensionalData
using Interpolations

@kwdef mutable struct Config 
    base_path::String
    path_to_recipes::String
    prefix_recipes::String
    target_dir::String
    diagnostics::Vector{String}
    variables::Vector{String}
    name_ref_period::String
    name_full_period::String
    name_obs_period::String
    models_project_name::String
    obs_data_name::String
    weights_variables::Union{Nothing, Dict{String, Dict{String,Number}}}
    weight_contributions::Union{Nothing, Dict{String, Number}}
end


function warnIfFlawedMetadata(attributes, filename)
    isFlawed = false;
    if "branch_time_in_parent" in keys(attributes) && isa(attributes["branch_time_in_parent"], String)
        @warn "Branch_time_in_parent is a string, excluded file:" filename
        isFlawed = true
    elseif "branch_time_in_child" in keys(attributes) && isa(attributes["branch_time_in_child"], String)
        @warn "Branch_time_in_child is a string, excluded file:" filename
        isFlawed = true
    end
    return isFlawed
end


function buildCMIP5EnsembleMember(realizations, initializations, physics, n)
    function concat(elem, prefix)
        if !(elem isa Vector)
            return [prefix * string(elem) for _ in range(1, n)]
        else 
            return map(x -> prefix * string(x), elem)
        end
    end
    rips = [concat(realizations, "r"), concat(initializations, "i"), concat(physics, "p")];
    variants = [join(k, "") for k in zip(rips...)]
    return variants
end


"""
    updateMetadata!(
        meta::Dict{String, Union{Array, String}}, 
        data::DimArray, 
        indices_ignored::Array
    )

Update metadata 'meta' s.t. data of ignored files is removed and attributes
that were only present in some files/models are set to missing. Further the key 
'full_model_names' is added which contains for every file/model the unique 
identifier consisting of variant_label, grid_label and model_name. 
"""
function updateMetadata!(
    meta::Dict{String, Union{Array, String, Dict}}, 
    data::DimArray
)
    model_dim = length(dims(data, :model))
    for key in keys(meta)
        values = meta[key]; 
        # add missing values for the last added files 
        n =  model_dim - length(values)
        for _ in range(1, n)
            push!(values, missing)
        end
        #println("i: " * string(i) * " " * string(length(meta[key])))
        # if none was missing and all have the same value, just use a string
        if !any(ismissing, values) && length(unique(values)) == 1
            meta[key] = string(values[1])
        end
    end
    # for model data only
    if !isempty(get(meta, "source_id", get(meta, "model_id", "")))
        key_model_name = getCMIPModelsKey(meta);
        meta["full_model_names"] = getUniqueModelIds(meta, key_model_name)
        # add mapping from model (ensemble) names to indices in metadata arrays
        meta["ensemble_names"] = Array(dims(data, :model))
        meta["indices_map"] = getIndicesMapping(meta["ensemble_names"]);
    end

end


"""
    appendValuesDicts(val1, val2)

Combine values of two different dictionaries iteratively. If both values are 
Vectors, they should be identical and only one of them is added. If they
aren't identical, a warning is triggered and both vectors are concatenated 
into one big Vector.

# Arguments:
- `val1`: a Vector, a Dictionary or a single value (e.g. String, Number, etc.)
- `val2`: a Vector, a Dictionary or a single value (e.g. String, Number, etc.)
"""
function appendValuesDicts(val1, val2)
    if isa(val1, Vector) && isa(val2, Vector) & !(isempty(val1) || isempty(val2))
        if val1 != val2 
            @warn "Two arrays merged that weren't identical! (usuallly in metadata)"
            @warn val1
            @warn val2
            return vcat(val1, val2)
        else
            return val1 
        end
    elseif isa(val1, Vector)
        return push!(val1, val2)
    elseif isa(val2, Vector)
        return push!(val2, val1)

    elseif isa(val1, Dict) && isa(val2, Dict)
        if val1 != val2
            @warn "Two different dictionaries (check indices_map)"
            @warn val1 
            @warn val2
            return vcat(val1, val2)
        else
            return val1
        end

    elseif isa(val1, Dict) || isa(val2, Dict)
        @warn val1 
        @warn val2
        throw(ArgumentError("Dictionary merged with something else!"))
    else
        return [val1, val2]
    end
end

"""
    buildPathsToVarData(config::Config, name_time_period::String)

Returns a mapping from diagnostic (e.g. CLIM) to variable (e.g. tas) to the
paths where the data is stored. The assumed structure of the data is this: 
config.base_path contains the path to the directory that contains a directory 
called 'preproc', which contains a directory for every time period considered 
(e.g. historical1, historical) which in turn contains a directory for every 
combination of variable and diagnostic (e.g., tas_CLIM).
"""
function buildPathsToVarData(config::Config, name_time_period::String)
    var_to_path = Dict{String, Dict{String, String}}();
    for diagnostic in config.diagnostics
        var_to_path[diagnostic] = Dict();
        for var in config.variables
            var_diagnostic = var * "_" * diagnostic
            var_to_path[diagnostic][var] = joinpath(
                config.base_path,
                "preproc",
                name_time_period,
                var_diagnostic
            )
        end
    end
    return var_to_path
end


function validateConfig(path_config::String)
    config_yaml = YAML.load_file(path_config);
    weights_variables = get(config_yaml, "weights_variables", nothing);
    if !isnothing(weights_variables)
        weights_variables = convert(
            Dict{String, Dict{String, Number}}, 
            weights_variables
        )
    end
    # one of base_path and path_to_recipes must be provided
    base_path = get(config_yaml, "base_path", "")
    path_to_recipes = get(config_yaml, "path_to_recipes", "")
    if isempty(base_path) && isempty(path_to_recipes)
        msg = "Either 'base_path' (path to directory that contains preproc-directory) or 'path_to_recipes' (path to directory that contains for every variable a subfolder which in turn contain the preproc-directories.)";
        throw(ArgumentError(msg))
    end
    config = Config(
        base_path = base_path,
        path_to_recipes = path_to_recipes,
        prefix_recipes = get(config_yaml, "prefix_recipes", ""),
        target_dir = joinpath(config_yaml["target_dir"], getCurrentTime()),
        
        diagnostics = config_yaml["diagnostics"],
        variables = config_yaml["variables"],
        
        name_ref_period = get(config_yaml, "name_ref_period", ""),
        name_full_period = get(config_yaml, "name_full_period", ""),
        name_obs_period = get(config_yaml, "name_obs_period", ""),
        
        models_project_name = config_yaml["models_project_name"],
        obs_data_name = get(config_yaml, "obs_data_name", ""),
        
        weights_variables = weights_variables,
        weight_contributions = get(config_yaml, "weight_contributions", nothing)
    ) 
    # TODO: add checks consistency, paths for all specified variables there and exist, etc.
    return config
end


"""
    keepMetadataSubset!(meta::Dict, indices::Vector{Number})

# Arguments:
- `meta::Dict`: metadata dictionary
- `indices::Vector{Number}`: indices of data to be kept
"""
function keepMetadataSubset!(meta::Dict, indices::Vector{Int64})
    attributes = filter(x -> meta[x] isa Vector, keys(meta));
    for key in attributes
        meta[key] = meta[key][indices];
    end
end

"""
    keepModelSubset(data::Dict{String, DimArray}, shared_models::Vector{String})

Retain data only from models in `shared_models`. Takes care of metadata.
"""
function keepModelSubset(data::DimArray, shared_models::Vector{String})
    indices = findall(m -> m in shared_models, data.metadata["full_model_names"]);
    @assert length(indices) == length(shared_models)
    keepMetadataSubset!(data.metadata, indices);
    data = data[model = indices];
    data.metadata["indices_map"] = getIndicesMapping(data.metadata["ensemble_names"])
    return data
end


"""
    updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)

Vectors in metadata 'meta' have to refer to different models. These are now summarized such that
each vector only contains N entries where N is the number of Ensembles/Models (i.e. without the unique ensemble members).
If the metadata for the ensemble members of a Model/Ensemble differ across the members, the respective 
entry in the vector will be a vector itself. 
"""
function updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)
    meta_new = filter(((k,v),) -> !(v isa Vector), meta);    
    attributes = filter(x -> meta[x] isa Vector, keys(meta));
    attribs_diff_across_members = [];
    for key in attributes
        for model in dims(grouped_data, :model)
            indices = meta["indices_map"][model]
            vals = get!(meta_new, key, [])       
            val_ensemble = unique(meta[key][indices]);
            if length(val_ensemble) != 1
                push!(vals, val_ensemble)
                push!(attribs_diff_across_members, key)
            else
                push!(vals, val_ensemble[1])
            end
        end
    end
    if !isempty(attribs_diff_across_members)
        @warn "metadata attributes that differ across ensemble members (ok for some!)" unique(attribs_diff_across_members)
        @assert !(getCMIPModelsKey(meta) in attribs_diff_across_members)
    end
    return meta_new
end


"""
    computeInterpolatedWeightedQuantiles

This implementation follows the one used by Brunner et al.
"""
function computeInterpolatedWeightedQuantiles(quantiles, vals, weights=nothing)
    if isnothing(weights)
        weights = ones(length(vals));
    end
    indicesSorted = Array(sortperm(vals)); # gives indices of smallest to largest data point
    weightsSorted = weights[indicesSorted];
    weightedQuantiles = cumsum(weightsSorted) - 0.5 * weightsSorted;
    weightedQuantiles = reshape(weightedQuantiles, length(weightedQuantiles), 1);
    weightedQuantiles = (weightedQuantiles .- minimum(weightedQuantiles)) ./ maximum(weightedQuantiles);
    
    interp_linear = Interpolations.linear_interpolation(
        vec(weightedQuantiles), 
        vals[indicesSorted],
        extrapolation_bc=Interpolations.Line()
    );
     
    return interp_linear(quantiles)
end


"""
    renameModelDimsFromMemberToEnsemble(data::DimArray, dim_names::Vector{String})

# Arguments:
- `data`
- `dim_names`: names of dimensions to be changed, e.g. 'model', 'model1', etc.
"""
function renameModelDimsFromMemberToEnsemble(data::DimArray, dim_names::Vector{String})
    for dim in dim_names
        unique_members = dims(data, Symbol(dim))
        ensembles = map(x -> split(x, MODEL_MEMBER_DELIM)[1], unique_members)
        data = set(data, Symbol(dim) => ensembles)
    end
    return data
end