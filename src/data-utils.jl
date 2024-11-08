import YAML
using DimensionalData
using Interpolations


@kwdef struct DataID
    key::String
    variable::String
    statistic::String
    alias::String
    exp::String
    timerange::String
end

# Overload the Base.show method to print key-value pairs of DataID instances
function Base.show(io::IO, x::DataID)
    for field in fieldnames(DataID)
        value = getfield(x, field)        
        print(io, "$field=$value  ")
    end
end



@kwdef struct Data
    base_path::String
    ids::Vector{DataID}=[]
    data::Dict{String, DimArray}=Dict()
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
    performance_distances::DimArray
    independence_distances::DimArray
    Di::DimArray # generalized distances each model wrt performance
    Sij::DimArray # generalized distances between pairs of models
    wP::DimArray # normalized
    wI::DimArray # normalized
    w::DimArray # normalized
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


function buildCMIP5EnsembleMember(
    realizations::Vector, initializations::Vector, physics::Vector
)
    n = length(realizations)
    if (n != length(initializations)) || (n != length(physics))
        msg = "inconsistent input for building up CMIP5EnsembleMembers!"
        throw(ArgumentError(msg))
    end
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
        source_names::Vector{String},
        isModelData::Bool
    )

Update metadata 'meta' s.t. data of ignored files is removed and attributes
that were only present in some files/models are set to missing. Further the key 
'full_model_names' is added which contains for every file/model the unique 
identifier consisting of variant_label, grid_label (for CMIP6) and model_name. 
"""
function updateMetadata!(
    meta::Dict{String, Any}, 
    source_names::Vector{Union{Missing, String}},
    isModelData::Bool
)
    indices = findall(x -> !ismissing(x), source_names)
    for key in keys(meta)
        values = meta[key][indices]
        meta[key] = values  
        # if none was missing and all have the same value, just use a string
        if !any(ismissing, values) && length(unique(values)) == 1
            meta[key] = string(values[1])
        end
    end
    if isModelData
        included_models = Array{String}(source_names[indices])
        meta["full_model_names"] = getUniqueModelIds(meta, included_models)
        # add mapping from model (ensemble) names to indices in metadata arrays
        meta["ensemble_names"] = included_models
        meta["ensemble_indices_map"] = getIndicesMapping(included_models)
    end
    return nothing
end

"""
    joinMetadata(meta1::Dict{String, Any}, meta2::Dict{String, Any})

Join the metadata of data to be joint (e.g. when loading data for tos for
CMIP5 and CMIP6 from different locations). If keys are present in 'meta1' or 
'meta2' but not the other, missing values are added.
"""
function joinMetadata(meta1::Dict{String, Any}, meta2::Dict{String, Any})
    meta = Dict{String, Any}()
    n1 = length(meta1["full_model_names"])
    n2 = length(meta2["full_model_names"])
    keys_meta1 = keys(meta1)
    keys_meta2 = keys(meta2)
    keys_shared = collect(intersect(keys_meta1, keys_meta2))
    keys_uniq_m1 = filter(x -> !(x in keys_shared), keys_meta1)
    keys_uniq_m2 = filter(x -> !(x in keys_shared), keys_meta2)

    for k in keys_shared
        if k == "ensemble_map_indices"
            meta["ensemble_indices_map"] = Dict()
            for (k,vals) in meta2["ensemble_indices_map"]
                meta["ensemble_indices_map"][k] = vals .+ n1
            end
            for (k, vals) in meta1["ensemble_indices_map"]
                meta["ensemble_indices_map"][k] = deepcopy(vals)
            end
        else
            v1 = meta1[k]
            v2 = meta2[k]
            if isa(v1, String)
                if isa(v2, String)
                    if v1 == v2
                        meta[k] = v1
                    else
                        meta[k] = vcat(repeat([v1], outer=n1), repeat([v2], outer=n2))
                    end
                else
                    meta[k] = vcat(repeat([v1], outer=n1), v2)
                end
            elseif isa(v1, Vector)
                if isa(v2, String)
                    meta[k] = vcat(v1, repeat([v2], outer=n2))
                elseif isa(v2, Vector)
                    meta[k] = vcat(v1, v2)
                end
            end
        end
    end

    for keys_uniq in [keys_uniq_m1, keys_uniq_m2]
        for k in keys_uniq
            v = get(meta1, k, nothing)
            in_meta1 = true
            if isnothing(v)
                in_meta1 = false
                v = meta2[k]
                n_added = n1
                n = n2
            else
                n_added = n2
                n = n1
            end
            v_added = repeat([missing], outer=n_added)
            # be sure to add vectors in correct order (because v may refer to value of meta1 or meta2!)
            if isa(v, String)
                if in_meta1
                    meta[k] = vcat(repeat([v], outer=n), v_added)
                else
                    meta[k] = vcat(v_added, repeat([v], outer=n))
                end
            else
                if in_meta1
                    meta[k] = vcat(v, v_added)
                else
                    meta[k] = vcat(v_added, v)
                end
            end
        end
    end
    return meta
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
    if isa(val1, Vector) && isa(val2, Vector) #&& (!(isempty(val1) || isempty(val2)))
        if collect(skipmissing(val1)) != collect(skipmissing(val2)) 
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
            @warn "Two different dictionaries (check ensemble_indices_map)"
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

    elseif isa(val1, String) && isa(val2, String) && val1 == val2
        return val1
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
            indices = meta["ensemble_indices_map"][model]
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
    buildDataIDsFromConfigs(config_path::String)

# Arguments:
- `config_path`: path to directory that contains one or more yaml config
files. For the assumed structure of the config files, see: TODO.
"""
function buildDataIDsFromConfigs(config_path::String)
    paths_to_configs = filter(
        x -> isfile(x) && endswith(x, ".yml"), 
        readdir(config_path, join=true)
    )
    ids::Vector{DataID} = []
    for path_config in paths_to_configs
        config = YAML.load_file(path_config);
        data_all = config["diagnostics"]
        aliases = keys(data_all)
    
        for alias in aliases
            data = data_all[alias]["variables"]     
            for (k,v) in data
                variable, statistic = split(k, "_")
                if typeof(v["exp"]) <: String
                    experiment = v["exp"]
                else 
                    experiment = join(v["exp"], "-")
                end
                timerange = replace(get(v, "timerange", "full"), "/" => "-")
                id = join([variable, statistic, alias], "_")
                dataID = DataID(
                    key=id, 
                    variable=variable,
                    statistic=statistic,
                    alias=alias,
                    exp=experiment, 
                    timerange=timerange)
                push!(ids, dataID)
            end
        end
    end
    return ids
end


function applyDataConstraints!(ids::Vector{DataID}, subset::Dict{String, Vector{String}})   
    if !isempty(get(subset, "variables", Vector{String}()))
        keepVar(id::DataID) = any(var -> id.variable == var, subset["variables"])
        filter!(keepVar, ids)
    end
    if !isempty(get(subset, "aliases", Vector{String}()))
        keepTasks(id::DataID) = any(name -> id.alias == name, subset["aliases"])
        filter!(keepTasks, ids)
    end

    if !isempty(get(subset, "statistics", Vector{String}()))
        keepStats(id::DataID) = any(stat -> id.statistic == stat, subset["statistics"])
        filter!(keepStats, ids)
    end

    if !isempty(get(subset, "timeranges", Vector{String}()))
        keepTimeRange(id::DataID) = any(tr -> id.timerange == tr, subset["timeranges"])
        filter!(keepTimeRange, ids)
    end
end


function alignIDsFilteredData!(data::Data)
    actual_data = keys(data.data)
    filter!(x -> x.key in actual_data, data.ids)
end
