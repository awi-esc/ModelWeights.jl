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
that were only present in some files/models are set to missing. Further keys
are added: 'full_model_names' contains for every file/model the unique 
identifier consisting of variant_label, grid_label (for CMIP6) and model_name, 
'model_names' contains the names of the included models without the ensemble 
member identifiers and 'model_to_member_indices' is a dictionary mapping from 
model names to the indices of the respective ensemble members.

Arguments:
- `meta`:
- `source_names`: 
- `isModelData`:
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
    included_data = Array{String}(source_names[indices])
    if isModelData
        meta["member_names"] = getUniqueModelIds(meta, included_data)
        # add mapping from model (ensemble) names to indices in metadata arrays
        meta["model_names"] = included_data
        meta["model_to_member_indices"] = getIndicesMapping(included_data)
    else 
        meta["source_names"] = included_data
    end
    return nothing
end

"""
    joinMetadata(meta1::Dict{String, Any}, meta2::Dict{String, Any}, isModelData::Bool)

Join the metadata of data to be joint (e.g. when loading data for tos for
CMIP5 and CMIP6 from different locations). If keys are present in 'meta1' or 
'meta2' but not the other, missing values are added.
"""
function joinMetadata(meta1::Dict{String, Any}, meta2::Dict{String, Any}, isModelData::Bool)
    meta = Dict{String, Any}()
    n1 = isModelData ? length(meta1["member_names"]) : length(meta1["source_names"])
    n2 = isModelData ? length(meta2["member_names"]) : length(meta2["source_names"])
    keys_meta1 = keys(meta1)
    keys_meta2 = keys(meta2)
    keys_shared = collect(intersect(keys_meta1, keys_meta2))
    keys_uniq_m1 = filter(x -> !(x in keys_shared), keys_meta1)
    keys_uniq_m2 = filter(x -> !(x in keys_shared), keys_meta2)

    for k in keys_shared
        if k == "model_to_member_indices"
            meta["model_to_member_indices"] = Dict()
            for (k,vals) in meta2["model_to_member_indices"]
                meta["model_to_member_indices"][k] = vals .+ n1
            end
            for (k, vals) in meta1["model_to_member_indices"]
                meta["model_to_member_indices"][k] = deepcopy(vals)
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
            @warn "Two different dictionaries (check model_to_member_indices)"
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
    meta["model_to_member_indices"] = getIndicesMapping(meta["member_names"])
    return nothing
end



"""
    updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)

Vectors in metadata 'meta' refer to different models (members). 
These are now summarized such that each vector only contains N entries where N
is the number of Ensembles (i.e. without the unique ensemble members).
If the metadata for the ensemble members of a Model/Ensemble differ across the
members, the respective entry in the vector will be a vector itself. 
"""
function updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)
    meta_new = filter(((k,v),) -> !(v isa Vector), meta)
    meta_new["model_to_member_indices"] = Dict{String, Number}()
    attributes = filter(x -> meta[x] isa Vector, keys(meta))
    attribs_diff_across_members = [];
    # iterate over attributes that are vectors, thus different for the different 
    # models or ensembles
    for key in attributes
        for (i, model) in enumerate(dims(grouped_data, :model))
            indices = meta["model_to_member_indices"][model]
            vals = get!(meta_new, key, [])       
            val_model = meta[key][indices]
            if length(unique(val_model)) != 1
                push!(vals, val_model)
                push!(attribs_diff_across_members, key)
            else
                # members of current ensemble all share the same value
                push!(vals, val_ensemble[1])
            end
            meta_new["model_to_member_indices"][model] = i
        end
    end
    if !isempty(attribs_diff_across_members)
        # TODO ignore those that are defenitely expected to differ
        @warn "metadata attributes that differ across model members (ok for some!)" unique(attribs_diff_across_members)
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
                    timerange=timerange
                )
                push!(ids, dataID)
            end
        end
    end
    return unique(ids)
end

"""
    applyDataConstraints!(ids::Vector{DataID}, subset::Dict{String, Vector{String}})

Subset model ids so that only those with properties specified in 'subset' remain.

# Arguments
- `ids`: Vector of DataID instances.
- `subset`: Mapping from fieldnames of 'DataID' struct to Vector specifiying the
properties of which at least one must be present for an id to be retained.
"""
function applyDataConstraints!(ids::Vector{DataID}, subset::Dict{String, Vector{String}})   
    for field in fieldnames(DataID)
        constraints = get(subset, string(field), Vector{String}()) # e.g. [historical, historical0]
        if !isempty(constraints)
            fn(id::DataID) = any(x -> getproperty(id, field) == x, constraints)
            filter!(fn, ids)
        end
    end
    return nothing
end


function indexData(data::Data, var_diagnostic_key::String)
    data_dict = filter(((k, v),) -> occursin(var_diagnostic_key, k), data.data)
    # check that only one dataset for combination of variable + diagnostic
    if length(data_dict) > 1
        var, diagnostic = split(var_diagnostic_key, "_")
        key = first(keys(data_dict))
        @warn "more than one dataset for $var and $diagnostic in model data, $key is taken!"
    end
    data = first(values(data_dict))
    return data
end


function computeDistancesAllDiagnostics(
    model_data::Data, obs_data::Data, var_diagnostic_keys::Vector{String}, forPerformance::Bool
)
    # compute performance/independence distances for all models and ensemble members
    diagnostics = unique(map(x -> split(x, "_")[2], var_diagnostic_keys))
    distances_all = []
    for diagnostic in diagnostics
        distances = []
        diagnostic_keys = filter(x -> endswith(x, "_" * diagnostic), var_diagnostic_keys)
        variables = map(x -> split(x, "_")[1], diagnostic_keys)
        for var in variables
            k = var * "_" * diagnostic
            models = indexData(model_data, k)
            
            if forPerformance
                observations = indexData(obs_data, k)
                if length(dims(observations, :source)) != 1
                    @warn "several observational datasets available for computing distances"
                end
                observations = observations[source=1]
                dists = getModelDataDist(models, observations)
            else
                dists = getModelDistances(models)
            end
            push!(distances, dists)
        end
        distances = cat(distances..., dims = Dim{:variable}(collect(variables)));
        push!(distances_all, distances)
    end
    return cat(distances_all..., dims = Dim{:diagnostic}(collect(diagnostics)));
end

function computeGeneralizedDistances(distances_all::DimArray, weights::DimArray, forPerformance::Bool)
    dimensions = forPerformance ? (:model,) : (:model1, :model2)
    norm = mapslices(Statistics.median, distances_all, dims=dimensions)
    normalized_distances =  DimArray(
        distances_all ./ norm, dims(distances_all), metadata = distances_all.metadata
    )
    distances = forPerformance ? 
        summarizeEnsembleMembersVector(normalized_distances, false) :
        averageEnsembleMatrix(normalized_distances, false);

    distances = mapslices(x -> x .* weights, distances, dims=(:variable, :diagnostic))
    return dropdims(
        sum(distances, dims=(:variable, :diagnostic)), dims=(:variable, :diagnostic)
    )
end

