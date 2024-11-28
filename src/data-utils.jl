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
    base_paths::Vector{String}
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
are added.
For model data:
    - 'member_names': vector that contains for every model a vector with the
    unique names of that model's members
    identifier consisting of variant_label, model_name and for CMIP6 models also grid_label.
    - 'model_names': vector whose length is the sum of the number of all models'
    members; it contains the model names for each unique model member, i.e. this
    vector will not unique if any model had several members
For observational data:
    - 'source_names': vector of data sources

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
        meta["member_names"] = getUniqueMemberIds(meta, included_data)
        meta["model_names"] = included_data
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
    n1 = isModelData ? length(vcat(meta1["member_names"]...)) : length(meta1["source_names"])
    n2 = isModelData ? length(vcat(meta2["member_names"]...)) : length(meta2["source_names"])
    keys_meta1 = keys(meta1)
    keys_meta2 = keys(meta2)
    keys_shared = collect(intersect(keys_meta1, keys_meta2))
    keys_uniq_m1 = filter(x -> !(x in keys_shared), keys_meta1)
    keys_uniq_m2 = filter(x -> !(x in keys_shared), keys_meta2)

    for k in keys_shared
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
            @warn "Two different dictionaries!"
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
    updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)

Vectors in metadata 'meta' refer to different models (members).
These are now summarized such that each vector only contains N entries where N
is the number of models (i.e. without the unique members).
If the metadata for members of a model differ across members, the respective
entry in the vector will be a vector itself.
"""
function updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)
    meta_new = filter(((k,v),) -> k=="member_names" || !(v isa Vector), meta)
    attributes = filter(x -> meta[x] isa Vector && x != "member_names", keys(meta))
    attribs_diff_across_members = [];
    # iterate over attributes that are vectors, thus different for the different
    # members or models
    for key in attributes
        for (i, model) in enumerate(dims(grouped_data, :model))
            indices = findall(x -> x==model, meta["model_names"])
            vals = get!(meta_new, key, [])
            val_model = meta[key][indices]
            if length(unique(val_model)) != 1
                push!(vals, val_model)
                push!(attribs_diff_across_members, key)
            else
                # members of current ensemble all share the same value
                push!(vals, val_model[1])
            end
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
    setLookupsFromMemberToModel(data::DimArray, dim_names::Vector{String})

Change the lookup values for the dimension 'member' to refer to the models, i.e.
they are not unique anymore. This is done in preparation to group the data by
the different models.

# Arguments:
- `data`: has at least dimensions in 'dim_names'
- `dim_names`: names of dimensions to be changed, e.g. 'member', 'member1'
"""
function setLookupsFromMemberToModel(data::DimArray, dim_names::Vector{String})
    n_dims = length(dim_names)
    for (i, dim) in enumerate(dim_names)
        unique_members = dims(data, Symbol(dim))
        models = map(x -> split(x, MODEL_MEMBER_DELIM)[1], unique_members)

        data = set(data, Symbol(dim) => models)
        new_dim_name = n_dims > 1 ? "model" * string(i) : "model"
        data = set(data, Symbol(dim) => Symbol(new_dim_name))
    end
    return data
end


"""
    buildDataIDsFromConfigs(config_paths::Vector{String})

# Arguments:
- `config_paths`: list of paths to directories that contain one or more yaml config
files. For the assumed structure of the config files, see: TODO.
"""
function buildDataIDsFromConfigs(config_paths::Vector{String})
    ids::Vector{DataID} = []
    for config_path in config_paths
        paths_to_configs = filter(
            x -> isfile(x) && endswith(x, ".yml"),
            readdir(config_path, join=true)
        )
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

                    # save the mapping between timeranges and aliases globally
                    TIMERANGE_TO_ALIAS[timerange] = alias
                    ALIAS_TO_TIMERANGE[alias] = timerange

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

    # check for compatibility of timerange and alias first
    timerange_constraints = get(subset, "timerange", Vector{String}())
    alias_constraints = get(subset, "alias", Vector{String}())
    timerangeOk(id::DataID) = any(x -> id.timerange == x, timerange_constraints)
    aliasOk(id::DataID) = any(x -> id.alias == x, alias_constraints)

    if !isempty(timerange_constraints) && !isempty(alias_constraints)
        filter!(x -> timerangeOk(x) || aliasOk(x), ids)
        if isempty(ids)
            msg = "Neither timeranges $(timerange_constraints) nor aliases $(alias_constraints) found in data!"
            throw(ArgumentError(msg))
        end
    elseif !isempty(timerange_constraints)
        filter!(timerangeOk, ids)
    elseif !isempty(alias_constraints)
        filter!(aliasOk, ids)
    end

    fields = filter(x -> !(x in [:key, :timerange, :alias]), fieldnames(DataID))
    for field in fields
        constraints = get(subset, string(field), Vector{String}()) # e.g. [historical, historical0]
        if !isempty(constraints)
            fn(id::DataID) = any(x -> getproperty(id, field) == x, constraints)
            filter!(fn, ids)
        end
    end
    return nothing
end


function indexData(data::Data, clim_var::String, diagnostic::String, ref_period::String)
    data_key = join([clim_var, diagnostic, ref_period], "_")
    return data.data[data_key]
end


function computeDistancesAllDiagnostics(
    model_data::Data,
    obs_data::Data,
    config::Dict{String, Number},
    ref_period::String,
    forPerformance::Bool
)
    # compute performance/independence distances for all model members
    var_diagnostic_keys = collect(keys(config))
    diagnostics = String.(unique(map(x -> split(x, "_")[2], var_diagnostic_keys)))
    distances_all = []
    for diagnostic in diagnostics
        distances = []
        diagnostic_keys = filter(x -> endswith(x, "_" * diagnostic), var_diagnostic_keys)
        variables = String.(map(x -> split(x, "_")[1], diagnostic_keys))
        for clim_var in variables
            models = indexData(model_data, clim_var, diagnostic, ref_period)

            if forPerformance
                observations = indexData(obs_data, clim_var, diagnostic, ref_period)
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
        distances = cat(distances..., dims = Dim{:variable}(String.(collect(variables))));
        push!(distances_all, distances)
    end
    return cat(distances_all..., dims = Dim{:diagnostic}(String.(collect(diagnostics))));
end


function computeGeneralizedDistances(distances_all::DimArray, weights::DimArray, forPerformance::Bool)
    dimensions = forPerformance ? (:member,) : (:member1, :member2)
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


function allcombinations(v...)
    combis = Vector{String}()
    for elems in Iterators.product(v...)
        push!(combis, join(elems, "_"))
    end
    return combis
end


"""
    verifyDataAndWeightInput(
    model_data::Data, obs_data::Data, weights::DimArray
)

Check that there is data for all keys of the provided (balance) weights and the given reference period 'ref_period'.

# Arguments:
- `data`:
- `keys_weights`:
- `ref_period`:
"""
function isValidDataAndWeightInput(data::Data, keys_weights::Vector{String}, ref_period::String)
    ids = map(x -> x.key, data.ids)
    keys_data = map(x -> x * "_" * ref_period, keys_weights)
    return all([k in ids for k in keys_data])
end


function getRefPeriodAsAlias(ref_period::String)
    return ref_period in keys(TIMERANGE_TO_ALIAS) ? TIMERANGE_TO_ALIAS[ref_period] : ref_period
end


function setRefPeriodInWeightsMetadata!(
    meta::Dict, ref_period_orig::String, ref_period_alias::String
)
    meta["ref_period_alias"] = ref_period_alias
    if ref_period_orig != ref_period_alias
        meta["ref_period_timerange"] = ref_period_orig
    else
        meta["ref_period_timerange"] = ALIAS_TO_TIMERANGE[ref_period_orig]
    end
    return nothing
end