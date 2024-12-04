import YAML
using DimensionalData
using Interpolations


@kwdef struct MetaAttrib
    variable::String
    statistic::String
    alias::String
    exp::String
    timerange::String
end


# Overload the Base.show method to print key-value pairs of MetaAttrib instances
function Base.show(io::IO, x::MetaAttrib)
    for field in fieldnames(MetaAttrib)
        value = getfield(x, field)
        print(io, "$field=$value ")
    end
end


@kwdef struct MetaData
    id::String
    attrib::MetaAttrib
    paths::Vector{String}
end


# Pretty print MetaData instances
function Base.show(io::IO, x::MetaData)
    println(io, "$(x.id) ($(x.attrib.timerange)) from:")
    for path in x.paths
        println(io, "\t$path")
    end
end


@kwdef struct Data
    meta::Vector{MetaData}=[]
    data::Dict{String, DimArray}=Dict()
end


"""
    joinDataObjects(data::Vector{Data})

# Return: A single Data-Object with the combined data from all Data-Objects 
in the input vector.
"""
function joinDataObjects(data::Vector{Data})
    meta_all = Vector{MetaData}()
    data_all = Dict{String, DimArray}()
    duplicate_ids = Vector{String}()
    for ds in data
        append!(meta_all, ds.meta)
        for k in keys(ds.data)
            if haskey(data_all, k) 
                push!(duplicate_ids, k)
                # TODO: handle duplicate ids (e.g, if for same metadata,
                # (different) data is loaded) 
                @warn "joint data with same id, one is overwritten!"
            end
            data_all[k] = copy(ds.data[k])
        end
    end
    return Data(data = data_all, meta = meta_all)
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

"""
    buildCMIP5EnsembleMember(
        realizations::Vector, initializations::Vector, physics::Vector
    )

Concatenate model settings to build ripf-abbreviations for CMIP5 models which 
do not have it in their metadata.

# Arguments
- `realizations`:
- `initializations`:
- `physics`:
"""
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
    getMetaAttribFromESMValToolConfigs(
        base_path_data::String,
        base_path_configs::String, 
        dir_per_var::Bool,
        subset::Dict{String, Vector{String}}
)

Read config files to determine which data shall be loaded, as defined by the 
variable, statistic, experiment and timerange. 

# Arguments:
- `config_paths`: list of paths to directories that contain one or more yaml config
files. For the assumed structure of the config files, see: TODO.
"""
function getMetaAttributesFromESMValToolConfigs(
    base_path_configs::String;
    subset::Union{Dict{String, Vector{String}}, Nothing} = nothing
)
    paths_configs = filter(
        x -> isfile(x) && endswith(x, ".yml"),
        readdir(base_path_configs, join=true)
    )
    meta_attrib = Vector{MetaAttrib}()
    for path_config in paths_configs
        config = YAML.load_file(path_config)
        data_all = config["diagnostics"]
        aliases = keys(data_all)

        for alias in aliases
            data = data_all[alias]["variables"]
            for (k,v) in data
                variable, statistic = String.(split(k, "_"))
                if typeof(v["exp"]) <: String
                    experiment = v["exp"]
                else
                    experiment = join(v["exp"], "-")
                end
                timerange = replace(get(v, "timerange", "full"), "/" => "-")
                meta = MetaAttrib(
                    variable = variable, 
                    statistic = statistic, 
                    alias = alias, 
                    exp = experiment,
                    timerange = timerange
                )
                if !(meta in meta_attrib) 
                    push!(meta_attrib, meta)
                end
            end
        end
    end
    if !isnothing(subset)
        applyDataConstraints!(meta_attrib, subset)
    end
    return meta_attrib
end



function getMetaDataFromYAML(
    path_config::String, 
    dir_per_var::Bool;
    subset::Union{Dict{String, Vector{String}}, Nothing} = nothing
)
    config = YAML.load_file(path_config)
    datasets = config["datasets"]
    base_path = get(config, "path_data", "")
    timerange_to_alias = config["timerange_to_alias"]
    meta_data = Vector{MetaData}()
    for ds in datasets
        path_data = joinpath(base_path, ds["base_dir"])
        subdir_constraints = get(ds, "subdirs", nothing) 
        for clim_var in ds["variables"]
            for stat in ds["statistics"]
                experiment = ds["exp"]                
                timeranges = get(ds, "timeranges", ["full"])
                for timerange in timeranges
                    alias = timerange == "full" ? experiment : 
                        get(timerange_to_alias, timerange, "unknown")
                    attrib = [MetaAttrib(
                        variable = clim_var, 
                        statistic = stat, 
                        alias = alias, 
                        exp = experiment,
                        timerange = timerange
                    )]
                    if !isnothing(subset)
                        applyDataConstraints!(attrib, subset)
                    end
                    if !isempty(attrib)
                        meta = buildMetaData(
                            attrib[1], path_data, dir_per_var; subdir_constraints
                        )     
                        append!(meta_data, meta)
                    end
                end
            end
        end
    end
    return meta_data
end




function buildMetaData(
    attrib::Union{MetaAttrib, Vector{MetaAttrib}},
    base_path_data::String,
    dir_per_var::Bool;
    subdir_constraints::Union{Vector{String}, Nothing} = nothing
)
    attributes = isa(attrib, MetaAttrib) ? [attrib] : attrib
    metadata = Vector{MetaData}()
    for attrib in attributes
        data_paths = buildPathsForMetaAttrib(
            base_path_data, attrib, dir_per_var; constraints=subdir_constraints
        )
        meta = MetaData(
            id = join([attrib.variable, attrib.statistic, attrib.alias], "_"),
            attrib = attrib,
            paths = data_paths
        )
        push!(metadata, meta)
    end
    return metadata
end

# if dir_per_var is true, directories at base_paths have subdirectories,
# one for each variable (they must contain '_VAR', e.g. '_tas'),
# otherwise base_paths are the paths to the directories that contain
# a subdirectory 'preproc'
function buildPathsForMetaAttrib(
    base_path::String, 
    attrib::MetaAttrib,
    dir_per_var::Bool;
    constraints::Union{Vector{String}, Nothing}=nothing
)
    base_paths = [base_path]
    if dir_per_var
        base_paths = filter(isdir, readdir(base_path, join=true))
        filter!(x -> occursin("_" * attrib.variable, x), base_paths)
    end
    if !isnothing(constraints) && !isempty(constraints)
        filter!(p -> any([occursin(name, p) for name in constraints]), base_paths)
    end
    data_paths = Vector{String}()
    for path in base_paths
        # Note: particular data structure assumed here!
        diagnostic = join([attrib.variable, attrib.statistic], "_")
        path_data = joinpath(path, "preproc", attrib.alias, diagnostic)
        if !isdir(path_data)
            @warn "No data found at:" path_data
        else
            push!(data_paths, path_data)
        end
    end
    return data_paths
end


"""
    applyDataConstraints!(
    metadata::Vector{MetaData}, subset::Dict{String, Vector{String}}
)

Subset model ids so that only those with properties specified in 'subset' remain.

# Arguments
- `metadata`: Vector of MetaData instances.
- `subset`: Mapping from fieldnames of 'MetaAttrib' struct to Vector specifiying the
properties of which at least one must be present for an id to be retained.
"""
function applyDataConstraints!(
    meta_attributes::Vector{MetaAttrib}, subset::Dict{String, Vector{String}}
)
    # check for compatibility of timerange and alias first
    timerange_constraints = get(subset, "timerange", Vector{String}())
    alias_constraints = get(subset, "alias", Vector{String}())
    timerangeOk(attrib::MetaAttrib) = any(x -> attrib.timerange == x, timerange_constraints)
    aliasOk(attrib::MetaAttrib) = any(x -> attrib.alias == x, alias_constraints)

    if !isempty(timerange_constraints) && !isempty(alias_constraints)
        filter!(x -> timerangeOk(x) || aliasOk(x), meta_attributes)
        if isempty(meta_attributes)
            msg = "Neither timeranges $(timerange_constraints) nor aliases $(alias_constraints) found in data!"
            throw(ArgumentError(msg))
        end
    elseif !isempty(timerange_constraints)
        filter!(timerangeOk, meta_attributes)
    elseif !isempty(alias_constraints)
        filter!(aliasOk, meta_attributes)
    end

    # constraints wrt remaining attributes
    fields = filter(x -> !(x in [:timerange, :alias]), fieldnames(MetaAttrib))
    for field in fields
        constraints = get(subset, string(field), Vector{String}())
        if !isempty(constraints)
            fn(attrib::MetaAttrib) = any(x -> getproperty(attrib, field) == x, constraints)
            filter!(fn, meta_attributes)
        end
    end
    return nothing
end


function indexData(data::Data, clim_var::String, diagnostic::String, ref_period::String)
    data_key = join([clim_var, diagnostic, ref_period], "_")
    return data.data[data_key]
end



"""
    computeDistancesAllDiagnostics(
        model_data::Data, 
        obs_data::Union{Nothing, Data}, 
        config::Dict{String, Number},
        ref_period::String,
        forPerformance::Bool    
    )

Compute RMSEs between models and observations or between predictions of models.

# Arguments:
- `model_data`:
- `obs_data`:
- `config`:
- `ref_period`:
- `for_performance`: true for distances between models and observations, false 
for distances between model predictions
"""
function computeDistancesAllDiagnostics(
    model_data::Data, 
    obs_data::Union{Nothing, Data}, 
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


"""
    computeGeneralizedDistances(
    distances_all::DimArray, weights::DimArray, for_performance::Bool
)

# Arguments:
- `distances_all`:
- `weights`:
- `for_performance`: 
"""
function computeGeneralizedDistances(
    distances_all::DimArray, weights::DimArray, for_performance::Bool
)
    dimensions = for_performance ? (:member,) : (:member1, :member2)
    norm = mapslices(Statistics.median, distances_all, dims=dimensions)
    normalized_distances =  DimArray(
        distances_all ./ norm, dims(distances_all), metadata = distances_all.metadata
    )
    distances = for_performance ? 
        summarizeEnsembleMembersVector(normalized_distances, false) :
        averageEnsembleMatrix(normalized_distances, false);

    distances = mapslices(x -> x .* weights, distances, dims=(:variable, :diagnostic))
    return dropdims(
        sum(distances, dims=(:variable, :diagnostic)), dims=(:variable, :diagnostic)
    )
end


"""
    allcombinations(v...)

Generate all possible combinations of input vectors, where each combination 
consists of one element from each input vector, concatenated as a string with
underscores separating the elements.

# Arguments
- `v...`: A variable number of input vectors.

# Returns
A vector of strings, where each string represents a unique combination of
elements from the input vectors, joined by underscores.

# Example
```jldoctest
julia> allcombinations(["tos", "tas"], ["CLIM"])


# output

2-element Vector{String}:
 "tos_CLIM"
 "tas_CLIM"
```
"""
function allcombinations(v...)
    combis = Vector{String}()
    for elems in Iterators.product(v...)
        push!(combis, join(elems, "_"))
    end
    return combis
end


"""
    isValidDataAndWeightInput(
        data::Data, keys_weights::Vector{String}, ref_period::String
    )

Check that there is data for all keys of the provided (balance) weights and 
the given reference period 'ref_period'.

# Arguments:
- `data`:
- `keys_weights`:
- `ref_period`:
"""
function isValidDataAndWeightInput(
    data::Data, keys_weights::Vector{String}, ref_period::String
)
    ids = map(x -> x.id, data.meta)
    keys_data = map(x -> x * "_" * ref_period, keys_weights)
    return all([k in ids for k in keys_data])
end

"""
    getTimerangeAsAlias(meta_attribs::Vector{MetaAttrib}, timerange::String)

Translate given timerange to corresponding alias in 'data_ids'.

# Arguments:
- `meta_attribs`:
- `timerange`:
"""
function getTimerangeAsAlias(meta_attribs::Vector{MetaAttrib}, timerange::String)
    attribs = filter(x -> x.timerange == timerange, meta_attribs)
    return isempty(attribs) ? nothing : attribs[1].alias
end


"""
    getAliasAsTimerange(meta_attribs::Vector{MetaAttrib}, alias::String)

Translate given alias to corresponding timerange in 'data_ids'.

# Arguments:
- `data_ids`:
- `alias`:
"""
function getAliasAsTimerange(meta_attribs::Vector{MetaAttrib}, alias::String)
    attribs = filter(x -> x.alias == alias, meta_attribs)
    return isempty(attribs) ? nothing : attribs[1].timerange
end


function getRefPeriodAsTimerangeAndAlias(
    meta_attribs::Vector{MetaAttrib}, ref_period::String
)
    alias = getTimerangeAsAlias(meta_attribs, ref_period) 
    # true : ref_period given as alias, false: ref_period given as timerange
    timerange = isnothing(alias) ? getAliasAsTimerange(meta_attribs, ref_period) : ref_period
    alias = isnothing(alias) ? ref_period : alias

    if (isnothing(alias) || isnothing(timerange))
        throw(ArgumentError("ref period $ref_period not given in data!"))
    end
    return (alias = alias, timerange = timerange)
end


function setRefPeriodInWeightsMetadata!(meta::Dict, alias::String, timerange::String)
    meta["ref_period_alias"] = alias
    meta["ref_period_timerange"] = timerange
    return nothing
end