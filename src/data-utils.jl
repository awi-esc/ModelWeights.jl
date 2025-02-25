import YAML
using DimensionalData
using Interpolations


@enum LEVEL MODEL=0 MEMBER=1

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
        if !isempty(value)
            print(io, "$field=$value ")
        end
    end
end


@kwdef struct MetaData
    id::String
    attrib::MetaAttrib
    paths::Vector{String}
end
# Pretty print MetaData instances
function Base.show(io::IO, x::MetaData)
    println(io, "::$(typeof(x)): $(x.id) (timerange: $(x.attrib.timerange), experiment: $(x.attrib.exp))")
end


@kwdef struct Data
    meta::MetaData
    data::DimArray
end
# Pretty print Data instances
function Base.show(io::IO, x::Data)
    println(io, "::$(typeof(x)): size: $(size(x.data)) id: $(x.meta.id)")
end

@kwdef struct ConfigWeights
    performance::Dict{String, Number}=Dict()
    independence::Dict{String, Number}=Dict()
    sigma_performance::Number=0.5
    sigma_independence::Number=0.5
    ref_perform_weights::String=""
    ref_indep_weights::String=""
    target_path::String=""
end


@kwdef struct Weights
    performance_distances::DimArray
    independence_distances::DimArray
    Di::DimArray # generalized distances each model wrt performance
    Sij::DimArray # generalized distances between pairs of models
    wP::DimArray # normalized
    wI::DimArray # normalized
    w::DimArray # normalized
    w_members::DimArray # weights distributed evenly across resp. model members
    config::ConfigWeights # metadata
end

# Pretty print Weights
function Base.show(io::IO, x::Weights)
    println(io, "::$(typeof(x)):")
    for m in dims(x.w, :model)
        println(io, "$m: $(round(x.w[model = At(m)], digits=3))")
    end
end



"""
    buildCMIP5EnsembleMember(
        realizations::Vector, initializations::Vector, physics::Vector
    )

Concatenate model settings to build ripf-abbreviations for CMIP5 models which 
do not have it in their metadata. Return a vector of strings.

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

Update metadata 'meta' s.t. attributes that were only present in some files/models 
are set to missing. Further key-value pairs are added concerning the data sources:
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
- `meta`: is the metadata dictionary of a DimArray
- `source_names`:
- `isModelData`:
"""
function updateMetadata!(
    meta_dict::Dict{String, Any},
    source_names::Vector{String},
    is_model_data::Bool
)
    sort_indices = sortperm(source_names)
    for key in keys(meta_dict)
        values = meta_dict[key][sort_indices]
        meta_dict[key] = values
        # if none was missing and all have the same value, just use a string
        if !any(ismissing, values) && length(unique(values)) == 1
            meta_dict[key] = string(values[1])
        end
    end
    # in some cases, the metadata model names do not match the model names as retrieved from the filenames 
    included_data = fixModelNamesMetadata(source_names[sort_indices])
    if is_model_data        
        member_ids = getUniqueMemberIds(meta_dict, included_data)
        meta_dict["member_names"] = member_ids
        meta_dict["model_names"] = included_data
        
        # if just data from one file is loaded, meta_dict["model_id"] is a string
        # (and we leave it as a string)
        get!(meta_dict, "model_id", "")
        if isa(meta_dict["model_id"], String)
            meta_dict["model_id"] = fixModelNamesMetadata([meta_dict["model_id"]])[1]
        else
            indices_non_missing = findall(map(x -> !ismissing(x), meta_dict["model_id"]))
            names = String.(meta_dict["model_id"][indices_non_missing])
            fixed_models = fixModelNamesMetadata(names)
            meta_dict["model_id"][indices_non_missing] = fixed_models
        end
        meta_dict["physics"] = getPhysicsFromMembers(member_ids)
    else
        meta_dict["source_names"] = included_data
    end
    return nothing
end


"""
    joinMetadata(
    meta1::Dict{String, Any}, meta2::Dict{String, Any}, isModelData::Bool
)

Join the metadata of data to be joint (e.g. when loading data for tos for
CMIP5 and CMIP6 from different locations). If keys are present in 'meta1' or 
'meta2' but not the other, missing values are added.

# Arguments:
- `meta1`:
- `meta2`:
- `is_model_data`: true for model data, false for observational data
"""
function joinMetadata(
    meta1::Dict{String, Any}, meta2::Dict{String, Any}, is_model_data::Bool
)
    meta = Dict{String, Any}()
    n1 = is_model_data ? length(meta1["member_names"]) : length(meta1["source_names"])
    n2 = is_model_data ? length(meta2["member_names"]) : length(meta2["source_names"])
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

# Arguments:
- `meta`:
- `grouped_data`:
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
        @debug "metadata attributes that differ across model members (ok for some!)" unique(attribs_diff_across_members)
    end
    return meta_new
end


"""
    computeInterpolatedWeightedQuantiles(
        quantiles::Vector{<:Number},
        vals::Vector;
        weights=nothing    
    )

This implementation follows the one used by Brunner et al.

# Arguments:
- `quantiles`: 
- `vals`:
- `weights`:
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

        data = DimensionalData.set(data, Symbol(dim) => models)
        new_dim_name = n_dims > 1 ? "model" * string(i) : "model"
        data = DimensionalData.set(data, Symbol(dim) => Symbol(new_dim_name))
    end
    return data
end


"""
    getMetaAttributesFromESMValToolConfigs(
        base_path_configs::String;
        constraint::Union{Dict, Nothing} = nothing
)

Read config files (ESMValTool recipes) to determine which data shall be loaded
in general (datapaths are not handled here), as defined by variable, statistic, 
experiment and timerange/alias.

# Arguments:
- `base_path_configs`: paths to directory that contain one or more yaml configs,
which may be ESMValTool recipes.
- `constraint`:

# Returns: A Vector of `MetaAttrib`-Objects.
"""
function getMetaAttributesFromESMValToolConfigs(
    base_path_configs::String;
    constraint::Union{Dict, Nothing} = nothing
)
    paths_configs = filter(
        x -> isfile(x) && endswith(x, ".yml"),
        readdir(base_path_configs, join=true)
    )
    meta_attribs = Vector{MetaAttrib}()
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
                if !(meta in meta_attribs) 
                    push!(meta_attribs, meta)
                end
            end
        end
    end
    if !isnothing(constraint)
        applyDataConstraints!(meta_attribs, constraint)
    end
    return meta_attribs
end


function addMetaData!(meta_dict::Dict{String, MetaData}, meta::MetaData)
    id = getMetaDataID(meta.attrib)
    if haskey(meta_dict, id)
        paths = meta_dict[id].paths
        append!(paths, meta.paths)
        meta_dict[id] = MetaData(
            id = id,
            attrib = meta.attrib,
            paths = unique(paths)
        )
    else
        meta_dict[id] = meta
    end
    return nothing
end


# for configuration with yaml file
function get_required_fields_config(ds::Dict)
    data = Dict(
        "base_dir" => get(ds, "base_dir", nothing),
        "exp" => get(ds, "exp", nothing),
        "variables" => get(ds, "variables", nothing),
        "statistics" => get(ds, "statistics", nothing)
    )
    if any(isnothing.(values(data)))
        # for surface altitude no statistics are computed, then we assume statistics="none"
        if isnothing(data["statistics"]) && data["variables"] == ["orog"]
            data["statistics"] = ["none"]
        else
            msg = "Config yaml file must specify values for the following required keys: 'exp' (experiment), 'base_dir' (path to data directory), 'variables' and 'statistics'."
            throw(ArgumentError(msg))
        end
    end
    return data
end


"""
    get_optional_fields_config(ds::Dict, timerange_aliases_dict::Dict)

Fill optional fields for a dataset in config file with default values. Given 
Returned timeranges and aliases have the same length and correspond to one another.

# Arguments:
- `ds`:
- `timerange_aliases_dict`:
"""
function get_optional_fields_config(ds::Dict, timerange_aliases_dict::Dict)
    aliases_timerange_vec = [(tr=k, alias=v) for (k,v) in timerange_aliases_dict]
    aliases_timerange_dict = Dict{String, String}()
    for elem in aliases_timerange_vec
        aliases_timerange_dict[elem.alias] = elem.tr
    end

    timeranges = unique(get(ds, "timeranges", Vector{String}()))
    full_included = "full" in timeranges
    filter!(x -> x != "full", timeranges)
    aliases =  unique(get(ds, "aliases", Vector{String}()))

    if !isempty(timeranges) && !isempty(aliases)
        # make sure that same data is not loaded once for timerange and once 
        # for corresponding alias
        tr_as_aliases = [get(timerange_aliases_dict, tr, nothing) for tr in timeranges]
        if any(isnothing.(tr_as_aliases))
            unknowns = timeranges[findall(x -> isnothing(x), tr_as_aliases)]
            throw(ArgumentError("Timeranges $unknowns aren't in timerange to alias dictionary in config file!"))
        end
        aliases_temp = filter(x -> !(x in tr_as_aliases), aliases)
        if !isempty(aliases_temp)
            aliases_as_tr = map(a -> filter(x -> x.alias==a, aliases_timerange_dict), aliases_temp)
            if any(x -> length(x) != 1, aliases_as_tr)
                throw(ArgumentError("Unknown alias according to timerange alias dictionary in config file!"))
            end
        else
            aliases_as_tr = Vector{String}()
        end
        aliases = vcat(tr_as_aliases, aliases_temp)
        timeranges = vcat(timeranges, aliases_as_tr)
    elseif isempty(timeranges)
        timeranges = [get(aliases_timerange_dict, a, nothing) for a in aliases]
        if any(isnothing.(timeranges))
            throw(ArgumentError("Unknown alias according to timerange alias dictionary in config file!"))
        end
    elseif isempty(aliases)
        aliases = [get(timerange_aliases_dict, tr, nothing) for tr in timeranges]
        if any(isnothing.(aliases))
            throw(ArgumentError("Unknown timerange according to timerange alias dictionary in config file!"))
        end
    end
    if  full_included || (isempty(timeranges) && isempty(aliases))
        push!(timeranges, "full")
        push!(aliases, ds["exp"])
    end
    return Dict(
        "timeranges" => timeranges,
        "aliases" => aliases,
        "models" => get(ds, "models", Vector{String}()),
        "projects" => get(ds, "projects", Vector{String}()),
        "subdirs" => get(ds, "subdirs", Vector{String}()),
        "dir_per_var" => get(ds, "dir_per_var", true)
    )
end

"""
If data is constraint by provided argument when loading, the argument takes
precedence over the given value inside the config yaml file.
"""
function setConstraintVal!(ds_config::Dict, cs_arg::Dict)
    for (field, val) in cs_arg
        ds_config[field] = val
    end
    return nothing
end

# function matchAliasTimerange!(constraint::Dict, timerange_to_alias::Dict)
#     if !isempty(aliases)
#         # use subset of aliases as given by input argument constraint
#         filter!(x -> get(timerange_to_alias, x, nothing) in constraint["aliases"], constraint["timeranges"])
#         if isempty(constraint["timeranges"])
#             msg = "Given aliases do not match with considered timeranges for dataset: $ds !"
#             throw(ArgumentError(msg))
#         end
#     end
# end


"""
    getMetaDataFromYAML(
        path_config::String,
        is_model_data::Bool
        subset::Union{Dict{String, Vector{String}}, Nothing} = nothing
    )

Load data as specified in config file located at `path_config`. For constraints
that are specified in the config file as well as in the `subset` argument, 
the values of the latter have precedence over the former. The constraints given
in the argument `subset` are applied to EVERY dataset specified in the config 
file.

# Arguments:
- `path_config`: path to config yaml file specifying meta attributes and pathes of data
- `is_model_data`: true for model data, false for observational data
- `subset`: TODO
"""
function getMetaDataFromYAML(
    path_config::String, 
    is_model_data::Bool;
    arg_constraint::Union{Dict, Nothing} = nothing
)
    config = YAML.load_file(path_config)
    datasets = config["datasets"]
    base_path = get(config, "path_data", "")
    timerange_to_alias = get(config, "timerange_to_alias", Dict{String, String}())


    meta_data = Dict{String, MetaData}()
    for ds in datasets
        # get data from config file
        req_fields = get_required_fields_config(ds)
        optional_fields = get_optional_fields_config(ds, timerange_to_alias)
        ds_constraint = merge(req_fields, optional_fields)
        # potentially update with constraint from argument
        if !isnothing(arg_constraint)
            setConstraintVal!(ds_constraint, arg_constraint)
        end
        if length(ds_constraint["timeranges"]) != length(ds_constraint["aliases"]) || isempty(ds_constraint["timeranges"])
            throw(ArgumentError("Timeranges and aliases must have the same length and match"))
        end
        # NOTE: if the second arg is an absolute path, joinpath will ignore the 
        # first arg and just use the second
        path_data = joinpath(base_path, req_fields["base_dir"])
        for clim_var in ds_constraint["variables"]
            for stats in ds_constraint["statistics"]
                for idx in 1:length(ds_constraint["timeranges"])
                    timerange = ds_constraint["timeranges"][idx]
                    alias =  ds_constraint["aliases"][idx]
                    attribs = [MetaAttrib(
                        variable = clim_var, 
                        statistic = stats, 
                        alias = alias, 
                        exp = ds_constraint["exp"],
                        timerange = timerange
                    )]
                    meta = buildMetaData(
                        attribs[1], path_data, ds_constraint["dir_per_var"], is_model_data; 
                        constraint = ds_constraint
                    )
                    if !isempty(meta.paths)
                        addMetaData!(meta_data, meta)
                    end
                end
            end
        end
    end
    return meta_data
end


"""
    buildPathsToDataFiles(
        path_data::String,
        is_model_data::Bool;
        model_constraints::Vector{String} = Vector{String}(),
        project_constraints::Vector{String} = Vector{String}()
    )

# Arguments:
- `path_data`:
- `is_model_data`: true for model data false for observational data
- `model_constraints`:
- `project_constraints`:

# Returns: Vector of Strings containing paths to data files in `path_data` 
that were not filtered out by `subset`.
"""
function buildPathsToDataFiles(
    path_data::String,
    is_model_data::Bool;
    model_constraints::Vector{String} = Vector{String}(),
    project_constraints::Vector{String} = Vector{String}()
)
    if !isdir(path_data)
        throw(ArgumentError(path_data * " does not exist!"))
    end
    if isempty(project_constraints)
        project_constraints = is_model_data ? ["CMIP"] : ["ERA5"]
    end
    ncFiles = filter(
        x -> isfile(x) && endswith(x, ".nc"),
        readdir(path_data; join=true)
    )
    @debug "$(map(println, ncFiles))"
    paths_to_files = Vector{String}()
    # constrain files that will be loaded
    for file in ncFiles
        keep = !isempty(project_constraints) ? 
            any([occursin(name, file) for name in project_constraints]) : true
        if !keep
            @debug "exclude $file because of projects subset"
            continue
        end
        keep = !isempty(model_constraints) ? 
            applyModelConstraints(file, model_constraints) : true
        if keep
            push!(paths_to_files, file)
        else 
            @debug "exclude $file because of model subset"
        end
    end
    return paths_to_files
end


function getMetaDataID(attrib::MetaAttrib)
    return join([attrib.variable, attrib.statistic, attrib.alias], "_")
end


"""
    buildMetaData(
        attrib::Union{MetaAttrib, Vector{MetaAttrib}},
        base_path_data::String,
        dir_per_var::Bool,
        is_model_data::Bool;
        constraint::Dict = Dict()
    )

Create data paths from assumed underlying data structure and `base_path_data` 
for every `MetaAttrib`-object in `attrib` and return a Vector of `MetaData`-
objects.

# Arguments:
- `attrib`:
- `base_path_data`:
- `dir_per_var`:
- `is_model_data`:
- `constraint`:

# Return: Vector of `MetaData`-objects.
"""
function buildMetaData(
    attrib::MetaAttrib,
    base_path_data::String,
    dir_per_var::Bool,
    is_model_data::Bool;
    constraint::Dict = Dict()
)
    subdir_constraints = isempty(constraint) ? nothing : get(constraint, "subdirs", Vector{String}())
    paths_data = buildPathsForMetaAttrib(
        base_path_data, attrib, dir_per_var; subdir_constraints
    )
    paths_to_files = Vector{String}()
    for path_data in paths_data
        paths = buildPathsToDataFiles(
            path_data, is_model_data; 
            model_constraints = get(constraint, "models", Vector{String}()),
            project_constraints = get(constraint, "projects", Vector{String}())
        )
        append!(paths_to_files, paths)
    end
    # for observational data, experiment doesn't make sense
    if !is_model_data
        attrib = MetaAttrib(
            variable = attrib.variable, 
            statistic = attrib.statistic,
            alias = attrib.alias,
            exp = "",
            timerange = attrib.timerange
        )
    end
    return MetaData(
        id = getMetaDataID(attrib),
        attrib = attrib,
        paths = paths_to_files
    )    
end

# if dir_per_var is true, directories at base_paths have subdirectories,
# one for each variable (they must contain '_VAR', e.g. '_tas'),
# otherwise base_paths are the paths to the directories that contain
# a subdirectory 'preproc'
"""
    buildPathsForMetaAttrib(
        base_path::String, 
        attrib::MetaAttrib,
        dir_per_var::Bool;
        subdir_constraints::Union{Vector{String}, Nothing}=nothing
    )

# Arguments:
- `base_path`: base directory of stored data specified in `attrib`.
- `attrib`: meta attributes of data.
- `dir_per_var`: true if data of each climate variable is stored in a seperate directory.
- `subdir_constraints`: if given, paths must contain ANY of the given elements. Existing paths that don't are ignored.
"""
function buildPathsForMetaAttrib(
    base_path::String, 
    attrib::MetaAttrib,
    dir_per_var::Bool;
    subdir_constraints::Union{Vector{String}, Nothing}=nothing
)
    base_paths = [base_path]
    if dir_per_var
        base_paths = filter(isdir, readdir(base_path, join=true))
        filter!(x -> occursin("_" * attrib.variable, x), base_paths)
    end
    if !isnothing(subdir_constraints) && !isempty(subdir_constraints)
        filter!(p -> any([occursin(name, p) for name in subdir_constraints]), base_paths)
    end
    data_paths = Vector{String}()
    for p in base_paths
        # Note: particular data structure assumed here!
        diagnostic = isempty(attrib.statistic) ? attrib.variable : 
            join([attrib.variable, attrib.statistic], "_")
        path_data = joinpath(p, "preproc", attrib.alias, diagnostic)
        if !isdir(path_data)
            @warn "$path_data is not an existing directory!"
        else
            push!(data_paths, path_data)
        end
    end
    return data_paths
end


"""
    applyDataConstraints!(meta_attributes::Vector{MetaAttrib}, constraint::Dict)

Subset model ids so that only those with properties specified in 'constraint' remain.

# Arguments
- `meta_attributes`: Vector of `MetaAttrib` instances.
- `constraint`: Mapping from fieldnames of 'MetaAttrib' struct to Vector 
specifiying the properties of which at least one must be present for an id to 
be retained.
"""
function applyDataConstraints!(
    meta_attributes::Vector{MetaAttrib}, constraint::Dict
)
    timerange_constraints = get(constraint, "timeranges", Vector{String}())
    alias_constraints = get(constraint, "aliases", Vector{String}())
    timerangeOk(attrib::MetaAttrib) = any(x -> attrib.timerange == x, timerange_constraints)
    aliasOk(attrib::MetaAttrib) = any(x -> attrib.alias == x, alias_constraints)

    # timerange and alias don't have to match, it's sufficient if either timerange
    # or alias match
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
    # constraints wrt variables and statistics
    stats_constraints = get(constraint, "statistics", Vector{String}())
    vars_constraints = get(constraint, "variables", Vector{String}())
    stats_ok(attrib::MetaAttrib) = isempty(stats_constraints) || any(x -> attrib.statistic == x, stats_constraints)
    vars_ok(attrib::MetaAttrib) = isempty(vars_constraints) || any(x -> attrib.variable == x, vars_constraints)
    filter!(stats_ok, meta_attributes)
    filter!(vars_ok, meta_attributes)
    return nothing
end


"""
    applyModelConstraints(file::String, model_constraints::Vector{String})

Returns true if constraints in `model_constraints` are fulfilled, i.e. if the 
given path to a model (`file`) contains any model from `model_constraints`,
false otherwise.

# Arguments:
- `file`: path to a file storing model data
- `model_constraints`: strings that may contain only model name, e.g. 'MPI-ESM-P', 
or model_name and member id, e.g. 'MPI-ESM-P#r1i1p2' or model name, member id and 
grid, e.g. 'MPI-ESM-P#r1i1p2_gn'
"""
function applyModelConstraints(file::String, model_constraints::Vector{String})
    keep_file = false
    for model in model_constraints
        keep_file = searchModelInPaths(model, [file])
        if keep_file
            break
        end
    end
    return keep_file
end


"""
    indexData(data::Data, clim_var::String, diagnostic::String, alias::String)
"""
function indexData(data::Vector{Data}, clim_var::String, diagnostic::String, alias::String)
    fn(data::Data) = 
        data.meta.attrib.variable == clim_var && 
        data.meta.attrib.statistic == diagnostic &&
        data.meta.attrib.alias == alias
    
    df = filter(fn, data)
    if length(df) == 0
        throw(ArgumentError("No data contained for $clim_var, $diagnostic, $alias !"))
    elseif length(df) > 1
        @warn "more than one dataset given for $clima_var, $diagnostic, $alias."
    end
    return df[1].data
end


"""
    computeDistancesAllDiagnostics(
        model_data::Vector{Data}, 
        obs_data::Union{Nothing, Data}, 
        config::Dict{String, Number},
        ref_period_alias::String,
        forPerformance::Bool    
    )

Compute RMSEs between models and observations or between predictions of models 
for all variables and diagnostics for which weights are specified in `config`.

# Arguments:
- `model_data`:
- `obs_data`:
- `config`: dictionary mapping from string of form VARIABLE_DIAGNOSTIC to respective weight.
- `ref_period_alias`:
- `for_performance`: true for distances between models and observations, 
false for distances between model predictions
"""
function computeDistancesAllDiagnostics(
    model_data::Vector{Data}, 
    obs_data::Union{Nothing, Vector{Data}}, 
    config::Dict{String, Number},
    ref_period_alias::String,
    for_performance::Bool
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
            models = indexData(model_data, clim_var, diagnostic, ref_period_alias)

            if for_performance
                observations = indexData(obs_data, clim_var, diagnostic, ref_period_alias)
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
        # TODO: recheck the metadata here! standard_name should be converted 
        # to a vector?!
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
    dimensions = for_performance ? 
        (hasdim(distances_all, :member) ? :member : :model) : 
        (:member1, :member2)
    norm = mapslices(Statistics.median, distances_all, dims=dimensions)
    normalized_distances =  DimArray(
        distances_all ./ norm, dims(distances_all), metadata = distances_all.metadata
    )
    if for_performance 
        distances = hasdim(normalized_distances, :model) ? 
            normalized_distances :
            summarizeEnsembleMembersVector(normalized_distances, false)
    else
        distances = averageEnsembleMembersMatrix(normalized_distances, false);
    end
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
julia> ModelWeights.allcombinations(["tos", "tas"], ["CLIM"])
2-element Vector{String}:
 "tos_CLIM"
 "tas_CLIM"
```
"""
function allcombinations(v...)
    if any(isempty.(v))
        @warn "At least one input vector is empty -> empty vector returned!"
    end
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

Check that there is data for all keys of the provided (balance) weights and the given reference period 'ref_period'.

# Arguments:
- `data`:
- `keys_weights`:
- `ref_period_alias`:
"""
function isValidDataAndWeightInput(
    data::Vector{Data}, keys_weights::Vector{String}, ref_period_alias::String
)
    actual_ids_data = map(x -> x.meta.id, data)
    required_keys_data = map(x -> x * "_" * ref_period_alias, keys_weights)
    return all([k in actual_ids_data for k in required_keys_data])
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

Translate given alias for every `MetaAttrib` in `meta_attribs`to corresponding timerange.

# Arguments:
- `meta_attribs`:
- `alias`:
"""
function getAliasAsTimerange(meta_attribs::Vector{MetaAttrib}, alias::String)
    attribs = filter(x -> x.alias == alias, meta_attribs)
    return isempty(attribs) ? nothing : attribs[1].timerange
end

"""
    getRefPeriodAsTimerangeAndAlias(
        meta_attribs::Vector{MetaAttrib}, ref_period::String
    )

Return timerange and alias corresponding to `ref_period` which can be either be 
a timerange or an alias.

# Arguments:
- `meta_attribs`:
- `ref_period`: either alias or timerange
"""
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

"""
    setRefPeriodInWeightsMetadata!(meta::Dict, alias::String, timerange::String)
"""
function setRefPeriodInWeightsMetadata!(meta::Dict, alias::String, timerange::String)
    meta["ref_period_alias"] = alias
    meta["ref_period_timerange"] = timerange
    return nothing
end


function fixModelNamesMetadata(names::Vector{String})
    model_names = copy(names)
    for m in keys(MODEL_NAME_FIXES)
        indices = findall(model_names .== m)
        if !isempty(indices)
            model_names[indices] .= MODEL_NAME_FIXES[m]
        end
    end
    return model_names
end

"""
    computeAreaWeights(longitudes::Vector{<:Number}, latitudes::Vector{<:Number}; mask::Union{DimArray, Nothing}=nothing)

Compute the approximated, normalized area weights for each lon,lat-position in 
`data` as the cosine of their latitudes. Return a DimensionalData.DimArray of
same size as the (lon x lat)-grid of the input data.

# Arguments:
- `data`: on lon, lat grid, with dimension 'lat' containing latitudes
- `mask`: optional 2D-mask (lonxlat), if given, the area weights are set to 0 where the mask is 1
"""
function computeAreaWeights(longitudes::Vector{<:Number}, latitudes::Vector{<:Number}; mask::Union{DimArray, Nothing}=nothing)
    # cosine of the latitudes as proxy for grid cell area
    areaWeights = cos.(deg2rad.(latitudes));
    areaWeights = DimArray(areaWeights, Dim{:lat}(latitudes));

    areaWeightMatrix = repeat(areaWeights', length(longitudes), 1);  
    if !isnothing(mask)
        areaWeightMatrix = ifelse.(mask .== 1, 0, areaWeightMatrix); 
    end
    areaWeightMatrix = areaWeightMatrix./sum(areaWeightMatrix)
    return DimArray(areaWeightMatrix, (Dim{:lon}(longitudes), Dim{:lat}(latitudes)))
end


function getMask(orog_data::DimArray, mask_out_land::Bool)
    ocean_mask = orog_data .== 0
    indices_missing = findall(x -> ismissing(x), ocean_mask);
    ocean_mask[indices_missing] .= false;
    ocean_mask = DimArray(Bool.(ocean_mask.data), dims(ocean_mask))
    mask = mask_out_land ? ocean_mask : ocean_mask .== false
    return mask
end

function getLandMask(orog_data::DimArray)
    return getMask(orog_data, false)    
end

function getOceanMask(orog_data::DimArray)
    return getMask(orog_data, true)    
end