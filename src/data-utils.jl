import YAML
using DimensionalData
using Interpolations


@kwdef struct DataID
    key::String
    exp::String
    statistic::String
    variable::String
    timerange::String
    task::String
end


@kwdef struct Data
    base_path::String
    ids::Vector{DataID}=[]
    models::Dict{String, DimArray}=Dict()
    obs::Dict{String, DimArray}=Dict()
    #data instead of models& obs
end



@kwdef struct ConfigWeights
    performance::Dict{String, Number}=Dict()
    independence::Dict{String, Number}=Dict()
    sigma_performance::Number=0.5
    sigma_independence::Number=0.5
    ref_period::String=""
    target_dir::String=""
end


@kwdef struct ClimwipWeights
    performance_all::DimArray
    independence_all::DimArray
    performance::DimArray #normalize
    independence::DimArray #normalize
    overall::DimArray #weights
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
        data::DimArray
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
    keepMetadataSubset!(meta::Dict, indices::Vector{Number})

# Arguments:
- `meta`: metadata dictionary
- `indices`: indices of data to be kept
"""
function keepMetadataSubset!(meta::Dict, indices::Vector{Int64})
    attributes = filter(x -> meta[x] isa Vector, keys(meta));
    for key in attributes
        meta[key] = meta[key][indices];
    end
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
    computeInterpolatedWeightedQuantiles(quantiles, vals; weights=nothing)

This implementation follows the one used by Brunner et al.
"""
function computeInterpolatedWeightedQuantiles(
    quantiles::Vector{<:Number},
    vals::Vector;
    weights=nothing
)
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


"""
    buildDataIDsFromConfigs(paths_to_config_dir::String)

# Arguments:
- `path_to_config_dir`: path to directory that contains one or more yaml config
files. For the assumed structure of the config files, see: TODO.
"""
function buildDataIDsFromConfigs(path_to_config_dir::String)
    paths_to_configs = filter(
        x -> isfile(x) && endswith(x, ".yml"), 
        readdir(path_to_config_dir, join=true)
    )
    ids::Vector{DataID} = []
    for path_config in paths_to_configs
        config = YAML.load_file(path_config);
        data_all = config["diagnostics"]
        names = keys(data_all)
    
        for name in names
            data = data_all[name]["variables"]     
            for (k,v) in data
                variable, statistic = split(k, "_")
                if typeof(v["exp"]) <: String
                    experiment = v["exp"]
                else 
                    experiment = join(v["exp"], "-")
                end
                timerange = replace(get(v, "timerange", "full"), "/" => "-")
                id = join([experiment, statistic, variable, timerange, name], "_", "#")
                dataID = DataID(id, experiment, statistic, variable, timerange, name)
                push!(ids, dataID)
            end
        end
    end
    return ids
end


function applyDataConstraints!(ids::Vector{DataID}, subset::Dict{String, Vector{String}})   
    if !isnothing(get(subset, "variables", nothing))
        keepVar(id::DataID) = any(var -> id.variable == var, subset["variables"])
        filter!(keepVar, ids)
    end
    if !isnothing(get(subset, "aliases", nothing))
        keepTasks(id::DataID) = any(name -> id.task == name, subset["aliases"])
        filter!(keepTasks, ids)
    end

    if !isnothing(get(subset, "statistics", nothing))
        keepStats(id::DataID) = any(stat -> id.statistic == stat, subset["statistics"])
        filter!(keepStats, ids)
    end

    if !isnothing(get(subset, "timeranges", nothing))
        keepTimeRange(id::DataID) = any(tr -> id.timerange == tr, subset["timeranges"])
        filter!(keepTimeRange, ids)
    end
end