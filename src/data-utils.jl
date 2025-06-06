"""
    joinDataMaps(v::DataMap...)
"""
function joinDataMaps(v::DataMap...)
    result = DataMap()
    for dm in v
        warnIfIdenticalKeys(result, dm)
        result = merge(result, dm)
    end
    return result
end


"""
    writeDataToDisk(data, target_path::String)

Save `data` as Julia obj if `target_path` has ending '.jld2', otherwise save as binary.
"""
function writeDataToDisk(data, target_path::String)
    if endswith(target_path, ".jld2")
        jldsave(target_path; data = data)
    else
        serialize(target_path, data)
    end
    @info "saved data to: $(target_path)"
    return nothing
end


"""
    readDataFromDisk(target_path::String; variable::String="")

Load data from `target_path`. If `target_path` ends with '.jld2', `variable` 
must be specified, otherwise  data is assumed to be binary.
"""
function readDataFromDisk(target_path::String; variable::String = "")
    if endswith(target_path, ".jld2")
        if isempty(variable)
            throw(ArgumentError("To load .jld2 data, specify argument variable!"))
        end
        f = jldopen(target_path, "r")
        data = f[variable]
        close(f)
    else
        data = deserialize(target_path)
    end
    return data
end


"""
    buildCMIP5EnsembleMember(
        realizations::Vector, initializations::Vector, physics::Vector
    )

Concatenate model settings to build ripf-abbreviations for CMIP5 models which 
do not have it in their metadata. Return a vector of strings.
"""
function buildCMIP5EnsembleMember(
    realizations::Vector,
    initializations::Vector,
    physics::Vector,
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
    rips = [concat(realizations, "r"), concat(initializations, "i"), concat(physics, "p")]
    variants = [join(k, "") for k in zip(rips...)]
    return variants
end


"""
    setLookupsFromMemberToModel(data::YAXArray, dim_names::Vector{String})

Change the lookup values for the dimension 'member' to refer to the models, i.e.
they are not unique anymore. This is done in preparation to group the data by
the different models.

# Arguments:
- `data::YAXArray`: has at least dimensions in `dim_names`.
- `dim_names::Vector{String}`: names of dimensions to be changed, e.g. 'member', 
'member1' (would be changed to 'model', 'model1').
"""
function setLookupsFromMemberToModel(data::YAXArray, dim_names::Vector{String})
    n_dims = length(dim_names)
    for (i, dim) in enumerate(dim_names)
        unique_members = dims(data, Symbol(dim))
        models = map(x -> String.(split(x, MODEL_MEMBER_DELIM)[1]), unique_members)

        data = DimensionalData.set(data, Symbol(dim) => models)
        new_dim_name = n_dims > 1 ? "model" * string(i) : "model"
        data = DimensionalData.set(data, Symbol(dim) => Symbol(new_dim_name))
    end
    return data
end


"""
    updateMetadata!(
        meta::Dict{String, Any},
        source_names::Vector{String},
        is_model_data::Bool
    )

Update metadata `meta` s.t. attributes that were only present in some 
files/models are set to missing. 

Further key-value pairs are added concerning the data sources:
- 'member_names' (for model data only): vector that contains for every model a 
vector with the unique names of that model's members identifier consisting of 
variant_label, model_name and for CMIP6 models also grid_label.
- 'model_names': vector whose length is the sum of the number of all models' 
members (or vector of filenames of observational data); it contains the model 
names for each unique model member, i.e. this vector will not be unique if any 
model had several members.

Arguments:
- `source_names::Vector{String}`: model names (or filenames of observational 
data)as retrieved from the original files the data was loaded from. 
"""
function updateMetadata!(
    meta_dict::Dict{String,Any},
    source_names::Vector{String},
    is_model_data::Bool,
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
    meta_dict["model_names"] = included_data

    if is_model_data
        member_ids = getUniqueMemberIds(meta_dict, included_data)
        meta_dict["member_names"] = member_ids

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
    end
    return nothing
end


"""
    updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)

Summarize vectors in `meta`, refering to different models (members), such that 
each vector only contains N entries (N=number of models (i.e. without unique 
members)).

If the metadata for members of a model differ across members, the respective
entry in the vector will be a vector itself.
"""
function updateGroupedDataMetadata(
    meta::Dict,
    grouped_data::DimensionalData.DimGroupByArray,
)
    meta_new = filter(((k, v),) -> k == "member_names" || !(v isa Vector), meta)
    attributes = filter(x -> meta[x] isa Vector && x != "member_names", keys(meta))
    attribs_diff_across_members = []
    # iterate over attributes that are vectors, thus different for the different
    # members or models
    for key in attributes
        for (i, model) in enumerate(dims(grouped_data, :model))
            indices = findall(x -> x == model, meta["model_names"])
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
        @debug "metadata attributes that differ across model members (ok for some!)" unique(
            attribs_diff_across_members,
        )
    end
    return meta_new
end


"""
    getMetaAttributesFromESMValToolConfigs(
        base_path_configs::String;
        constraint::Union{Dict, Nothing} = nothing
)

Read variable, statistic, experiment and timerange/alias values from ESMValTool 
recipes stored at `base_path_configs` into a vector of Dictionaries storing the 
respective readoff values.

If constraint is given, a combination of 'statistic', 'experiment', etc. is 
only retained if each of the respective readoff values aligns with the
the specified constraints if provided.
"""
function getMetaAttributesFromESMValToolConfigs(
    base_path_configs::String;
    constraint::Union{Dict,Nothing} = nothing,
)
    paths_configs = filter(
        x -> isfile(x) && endswith(x, ".yml"),
        readdir(base_path_configs, join = true),
    )
    meta_attribs = Vector{Dict{String,Any}}()
    for path_config in paths_configs
        config = YAML.load_file(path_config)
        data_all = config["diagnostics"]
        aliases = keys(data_all)

        for alias in aliases
            data = data_all[alias]["variables"]
            for (k, v) in data
                variable, statistic = String.(split(k, "_"))
                if typeof(v["exp"]) <: String
                    experiment = v["exp"]
                else
                    experiment = join(v["exp"], "-")
                end
                timerange = replace(get(v, "timerange", "full"), "/" => "-")
                meta = createGlobalMetaDataDict(
                    variable,
                    experiment,
                    statistic,
                    alias,
                    timerange,
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


# for configuration with yaml file
function get_required_fields_config(ds::Dict)
    data = Dict(
        "base_dir" => get(ds, "base_dir", nothing),
        "exp" => get(ds, "exp", nothing),
        "variables" => get(ds, "variables", nothing),
    )
    if any(isnothing.(values(data)))
        msg = "Config yaml file must specify values for the following required keys: $(keys(data))."
        throw(ArgumentError(msg))
    end
    # for fixed variables, no statistics are computed!
    # may add more than 'orog' for which no warning is thrown.
    stats = get(ds, "statistics", nothing)
    if isnothing(stats)
        data["statistics"] = ["none"]
        if data["variables"] != ["orog"]
            @warn "No statistics specified for $(data["variables"])."
        end
    else
        data["statistics"] = stats
    end
    return data
end


"""
    get_optional_fields_config(ds::Dict, timerange_aliases_dict::Dict)

Fill optional fields for a dataset in config file with default values. 
Returned timeranges and aliases have the same length and correspond to one another.

# Arguments:
- `ds`:
- `timerange_aliases_dict`:
"""
function get_optional_fields_config(ds::Dict, timerange_aliases_dict::Dict)
    aliases_timerange_vec = [(tr = k, alias = v) for (k, v) in timerange_aliases_dict]
    aliases_timerange_dict = Dict{String,String}()
    for elem in aliases_timerange_vec
        aliases_timerange_dict[elem.alias] = elem.tr
    end

    timeranges = unique(get(ds, "timeranges", Vector{String}()))
    full_included = "full" in timeranges
    filter!(x -> x != "full", timeranges)
    aliases = unique(get(ds, "aliases", Vector{String}()))

    if !isempty(timeranges) && !isempty(aliases)
        # make sure that same data is not loaded once for timerange and once 
        # for corresponding alias
        tr_as_aliases = [get(timerange_aliases_dict, tr, nothing) for tr in timeranges]
        if any(isnothing.(tr_as_aliases))
            unknowns = timeranges[findall(x -> isnothing(x), tr_as_aliases)]
            throw(
                ArgumentError(
                    "Timeranges $unknowns aren't in timerange to alias dictionary in config file!",
                ),
            )
        end
        aliases_temp = filter(x -> !(x in tr_as_aliases), aliases)
        if !isempty(aliases_temp)
            aliases_as_tr =
                map(a -> filter(x -> x.alias == a, aliases_timerange_dict), aliases_temp)
            if any(x -> length(x) != 1, aliases_as_tr)
                throw(
                    ArgumentError(
                        "Unknown alias according to timerange alias dictionary in config file!",
                    ),
                )
            end
        else
            aliases_as_tr = Vector{String}()
        end
        aliases = vcat(tr_as_aliases, aliases_temp)
        timeranges = vcat(timeranges, aliases_as_tr)
    elseif isempty(timeranges)
        timeranges = [get(aliases_timerange_dict, a, nothing) for a in aliases]
        if any(isnothing.(timeranges))
            unknowns = aliases[findall(x -> isnothing(x), timeranges)]
            throw(
                ArgumentError(
                    "Unknown alias according to timerange alias dictionary in config file: $unknowns",
                ),
            )
        end
    elseif isempty(aliases)
        aliases = [get(timerange_aliases_dict, tr, nothing) for tr in timeranges]
        if any(isnothing.(aliases))
            throw(
                ArgumentError(
                    "Unknown timerange according to timerange alias dictionary in config file!",
                ),
            )
        end
    end
    if full_included || (isempty(timeranges) && isempty(aliases))
        # default is "full"
        push!(timeranges, "full")
        push!(aliases, ds["exp"])
    end
    # sanity check
    if (length(timeranges) != length(aliases)) || isempty(timeranges)
        msg = "Merging timeranges and aliases failed! After merging, they must have the same length and match!"
        throw(ArgumentError(msg))
    end

    # subset this dataset only 
    subset_level = lowercase(get(ds, "subset_shared", ""))
    subset_level =
        subset_level == "model" ? MODEL : (subset_level == "member" ? MEMBER : nothing)

    return Dict(
        "timeranges" => timeranges,
        "aliases" => aliases,
        "models" => get(ds, "models", Vector{String}()),
        "projects" => get(ds, "projects", Vector{String}()),
        "subdirs" => get(ds, "subdirs", Vector{String}()),
        "dir_per_var" => get(ds, "dir_per_var", true),
        "subset_shared" => subset_level,
    )
end


# """
# If data is constraint by provided argument when loading, the argument takes
# precedence over the given value inside the config yaml file.
# """
# function setConstraintVal!(ds_config::Dict, cs_arg::Dict)
#     for (field, val) in cs_arg
#         ds_config[field] = val
#     end
#     return nothing
# end


"""
    getMetaDataFromYAML(
        content::Dict, is_model_data::Bool; arg_constraint::Union{Dict, Nothing} = nothing
    )

Return metadata of data specified in `content` possibly constrained by values in `arg_constraint`.

For constraints that are specified in `content` as well as in the `arg_constraint` argument, 
the values of the latter have precedence over the former. The constraints given
in the argument `arg_constraint` are applied to EVERY dataset specified in the config 
file.

# Arguments:
- `content`: content of config yaml file specifying meta attributes and paths of data
- `is_model_data`: true for model data, false for observational data
- `arg_constraint`: TODO
"""
function getMetaDataFromYAML(
    content::Dict, is_model_data::Bool; arg_constraint::Union{Dict, Nothing} = nothing
)
    datasets = content["datasets"]
    base_path = get(content, "path_data", "")
    timerange_to_alias = get(content, "timerange_to_alias", Dict{String,String}())

    meta_data = Dict{String, Dict{String, Any}}()
    for ds in datasets
        meta_ds = Dict{String, Dict{String, Any}}()
        # get data from config file
        req_fields = get_required_fields_config(ds)
        optional_fields = get_optional_fields_config(ds, timerange_to_alias)
        ds_constraint = merge(req_fields, optional_fields)
        considered_keys = keys(ds_constraint)
        # potentially update with constraint from argument (has precedence over value in yaml file)
        if !isnothing(arg_constraint)
            #setConstraintVal!(ds_constraint, arg_constraint)
            for (field, val) in arg_constraint
                # if subset_shared is given as argument, it will subset 
                # considering ALL loaded datasets, not individual ones!
                if field != "subset_shared"
                    if !(field in considered_keys)
                        @warn "$field is not a valid key to subset data!"
                    else
                        ds_constraint[field] = val
                    end
                end
            end
        end

        # NOTE: if the second arg is an absolute path, joinpath will ignore the 
        # first arg and just use the second
        path_data = joinpath(base_path, req_fields["base_dir"])
        for clim_var in ds_constraint["variables"]
            for stats in ds_constraint["statistics"]
                for idx = 1:length(ds_constraint["timeranges"])
                    timerange = ds_constraint["timeranges"][idx]
                    alias = ds_constraint["aliases"][idx]

                    meta = createGlobalMetaDataDict(
                        clim_var,
                        ds_constraint["exp"],
                        stats,
                        alias,
                        timerange
                    )
                    # for observational data, experiment doesn't make sense
                    if !is_model_data
                        meta["_experiment"] = ""
                    end
                    meta["_paths"] = getPathsToData(
                        meta,
                        path_data,
                        ds_constraint["dir_per_var"],
                        is_model_data;
                        constraint = ds_constraint,
                    )
                    if !isempty(meta["_paths"])
                        meta_ds[meta["_id"]] = meta
                    else
                        @warn "No data found for $(meta["_id"])"
                    end
                end
            end
        end
        # Apply subset_shared if provided inside yaml file for individual datasets
        if !isnothing(ds_constraint) && !isnothing(get(ds_constraint, "subset_shared", nothing))
            filterPathsSharedModels!(meta_ds, ds_constraint["subset_shared"])
        end
        # merge new dataset with already loaded
        for (id, meta) in meta_ds
            if haskey(meta_data, id)
                meta_data[id]["_paths"] = mergeMetaDataPaths(meta_data[id], meta)
            else
                meta_data[id] = meta
            end
        end
    end
    # filter across all datasets
    if !isnothing(arg_constraint) && !isnothing(get(arg_constraint, "subset_shared", nothing))
        filterPathsSharedModels!(meta_data, arg_constraint["subset_shared"])
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

Build vector of strings containing paths to data files in `path_data` 
that were not filtered out by `model_constraints` or `project_constraints`.

# Arguments:
- `model_constraints::Vector{String}`:
- `project_constraints::Vector{String}`:
"""
function buildPathsToDataFiles(
    path_data::String,
    is_model_data::Bool;
    model_constraints::Vector{String} = Vector{String}(),
    project_constraints::Vector{String} = Vector{String}(),
)
    if !isdir(path_data)
        throw(ArgumentError(path_data * " does not exist!"))
    end
    # if isempty(project_constraints)
    #     project_constraints = is_model_data ? ["CMIP"] : ["ERA5"]
    # end
    ncFiles = filter(x -> isfile(x) && endswith(x, ".nc"), readdir(path_data; join = true))
    # constrain files that will be loaded
    function doIncludeFile(file)
        keep = !isempty(project_constraints) ?
            any([occursin(name, file) for name in project_constraints]) : true
        if keep
            keep = !isempty(model_constraints) ? 
                applyModelConstraints(file, model_constraints) : true
        end
        return keep
    end
    mask = map(f -> doIncludeFile(f), ncFiles)
    return ncFiles[mask]
end


"""
    getPathsToData(
        attribs::Dict{String, Any},
        base_path_data::String,
        dir_per_var::Bool,
        is_model_data::Bool;
        constraint::Union{Dict, Nothing}=nothing
    )

Return the paths to the data files in `base_path_data` taking into account `constraint`.
"""
function getPathsToData(
    attribs::Dict{String, <:Any},
    base_path_data::String,
    dir_per_var::Bool,
    is_model_data::Bool;
    constraint::Union{Dict, Nothing} = nothing,
)
    has_constraint = !isnothing(constraint) && !isempty(constraint)
    subdir_constraints = !has_constraint ? nothing : get(constraint, "subdirs", Vector{String}())
    paths_data = buildPathsForMetaAttrib(
        base_path_data, attribs, dir_per_var; subdir_constraints
    )
    paths_to_files = Vector{String}()
    for path_data in paths_data
        paths = has_constraint ?
            buildPathsToDataFiles(
                path_data,
                is_model_data;
                model_constraints = get(constraint, "models", Vector{String}()),
                project_constraints = get(constraint, "projects", Vector{String}()),
            ) : 
            buildPathsToDataFiles(path_data, is_model_data)
        append!(paths_to_files, paths)
    end
    return paths_to_files
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

Return paths to data specified by `attribs`. 

Assumed data structure: BASE/preproc/ALIAS/VAR_STAT where BASE=`base_path`, 
ALIAS=`attribs["_alias"]`,VAR=attribs["_variable"] and STAT=attribs["_statistic"]. 
If `dir_per_var` is true, BASE is path to any subdirectory of `base_path` whose 
name contains _VAR. If `subdir_constraints` is given, the subdirectory's name must further 
contain at least one entry in `subdir_constraints`.

# Arguments:
- `base_path`: base directory of stored data specified in `attrib`.
- `attribs`: meta attributes of data. Must have keys: '_variable', '_statistic', '_alias'.
- `dir_per_var`: true if data of each climate variable is stored in a seperate directory.
- `subdir_constraints`: if given, paths must contain ANY of the given elements. Existing paths that don't are ignored.
"""
function buildPathsForMetaAttrib(
    base_path::String,
    attribs::Dict{String, <:Any},
    dir_per_var::Bool;
    subdir_constraints::Union{Vector{String},Nothing} = nothing,
)
    base_paths = [base_path]
    if dir_per_var
        base_paths = filter(isdir, readdir(base_path, join = true))
        filter!(x -> occursin("_" * attribs["_variable"], x), base_paths)
    end
    if !isnothing(subdir_constraints) && !isempty(subdir_constraints)
        filter!(p -> any([occursin(name, p) for name in subdir_constraints]), base_paths)
    end
    data_paths = Vector{String}()
    for p in base_paths
        # NOTE: particular data structure assumed here!
        diagnostic =
            isempty(attribs["_statistic"]) ? attribs["_variable"] :
            join([attribs["_variable"], attribs["_statistic"]], "_")
        path_data = joinpath(p, "preproc", attribs["_alias"], diagnostic)
        if isdir(path_data)
            push!(data_paths, path_data)
        end
        # else @debug "$path_data is not an existing directory!"
    end
    return data_paths
end


"""
    applyDataConstraints!(meta_attributes::Vector{Dict{String, Any}}, constraint::Dict)

Subset entries in `meta_attributes` so that only those with properties specified 
in `constraint` remain.

# Arguments
- `constraint::Dict`: Mapping to vector specifiying the properties of which at 
least one must be present for an id to be retained.
"""
function applyDataConstraints!(meta_attributes::Vector{Dict{String,Any}}, constraint::Dict)
    timerange_constraints = get(constraint, "timeranges", Vector{String}())
    alias_constraints = get(constraint, "aliases", Vector{String}())
    timerangeOk(attrib::Dict{String,Any}) =
        any(x -> attrib["_timerange"] == x, timerange_constraints)
    aliasOk(attrib::Dict{String,Any}) = any(x -> attrib["_alias"] == x, alias_constraints)

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
    stats_ok(attrib::Dict{String,Any}) =
        isempty(stats_constraints) || any(x -> attrib["_statistic"] == x, stats_constraints)
    vars_ok(attrib::Dict{String,Any}) =
        isempty(vars_constraints) || any(x -> attrib["_variable"] == x, vars_constraints)
    filter!(stats_ok, meta_attributes)
    filter!(vars_ok, meta_attributes)
    return nothing
end


"""
    applyModelConstraints(
        path_model_data::String, model_constraints::Vector{String}
    )

Return true if constraints in `model_constraints` are fulfilled, i.e. if `path_model_data` 
contains any model from `model_constraints`, false otherwise.

`path_model_data` must follow the standard CMIP filename structure, 
i.e. <variable_id>_<table_id>_<source_id>_<experiment_id >_<member_id>_<grid_label>[_<time_range>].nc 
and `model_constraints` can contain models on level of member (e.g. 'MPI-ESM-P#r1i1p2', 
'MPI-ESM-P#r1i1p2_gn') or on level of model (e.g., 'MPI-ESM-P').


# Arguments:
- `path_model_data`: 
- `model_constraints`: strings that may contain only model name, e.g. 'MPI-ESM-P', 
or model_name and member id, e.g. 'MPI-ESM-P#r1i1p2' or model name, member id and 
grid, e.g. 'MPI-ESM-P#r1i1p2_gn'.
"""
function applyModelConstraints(path_model_data::String, model_constraints::Vector{String})    
    keep_file = false
    for model in model_constraints
        keep_file = searchModelInPaths(model, [path_model_data])
        if keep_file
            break
        end
    end
    return keep_file
end


function avgObsDatasets(observations::YAXArray)
    if !hasdim(observations, :model)
        throw(ArgumentError("Obs data requires dim :model when averaging across datasets."))
    end
    if length(dims(observations, :model)) > 1
        return mean(observations, dims=:model)[model=At("combined")]
    else
        return observations[model=1]
    end
end


function distancesData(model_data::DataMap, obs_data::DataMap, config::Dict{String, Number})
    return distancesData(model_data, obs_data,  activeDiagnostics(config))
end


"""
    distancesData(model_data::DataMap, obs_data::DataMap, diagnostics_ids::Vector{String})

Compute RMSEs between models and observations for `diagnostics`.

# Arguments:
- `diagnostics_ids::Vector{String}`:
"""
function distancesData(model_data::DataMap, obs_data::DataMap, diagnostics::Vector{String})
    ensureDiagnosticsAvailable(model_data, diagnostics, "MODEL")
    ensureDiagnosticsAvailable(obs_data, diagnostics, "OBSERVATIONAL")
    distances_all = Vector{YAXArray}(undef, length(diagnostics))
    for (i, key) in enumerate(diagnostics)
        models = model_data[key]
        observations = obs_data[key]
        distances_all[i] =  distancesData(models, observations)
        # TODO: recheck the metadata here! standard_name should be converted to a vector?!
    end
    return cat(distances_all..., dims = Dim{:diagnostic}(diagnostics))
end


"""
    distancesData(models::YAXArray, observations::YAXArray)

Compute the distance as the area-weighted RMSE between model predictions and observations.

If several observational datasets are present, the average across all is taken.
"""
function distancesData(models::YAXArray, observations::YAXArray)
    obs = hasdim(observations, :model) ? avgObsDatasets(observations) : observations
    # Make sure to use a copy of the data, otherwise, it will be modified by applying the mask!!
    # I think this is not necessary actually (the deepcopy) 
    # Write test to be sure!!
    # models = deepcopy(models)
    n = length(dims(models, :member))
    member_names = Vector{String}(undef, n)
    distances = Vector{Any}(undef, n)
    for (i, model_i) in enumerate(eachslice(models; dims = :member))
        name = dims(models, :member)[i]
        maskNbMissing = (ismissing.(obs) + ismissing.(model_i)) .> 0 # observations or model is missing (or both)
        maskedObs = ifelse.(maskNbMissing .> 0, 0, obs)

        maskedModel = ifelse.(maskNbMissing .> 0, 0, model_i)
        mse = areaWeightedRMSE(maskedModel, maskedObs, maskNbMissing)

        member_names[i] = name
        distances[i] = mse
    end
    return YAXArray((Dim{:member}(member_names),), distances, models.properties)
end


function distancesModels(model_data::DataMap, config::Dict{String, Number})
    return distancesModels(model_data, activeDiagnostics(config))
end


function distancesModels(model_data::DataMap, diagnostics::Vector{String})
    ensureDiagnosticsAvailable(model_data, diagnostics, "MODEL")
    distances_all = Vector{YAXArray}(undef, length(diagnostics))
    for (i, key) in enumerate(diagnostics)
        models = model_data[key]
        distances_all[i] = distancesModels(models)
    end
    return cat(distances_all..., dims = Dim{:diagnostic}(diagnostics))
end


"""
    distancesModels(modelData::YAXArray)

Compute the area weighted RMSE between model predictions for each pair of models.

# Arguments:
- `modelData::YAXArray`: must have dimensions 'lon', 'lat', 'model' or 'member'.
"""
function distancesModels(model_data::YAXArray)
    # Make sure to use a copy of the data, otherwise, it will be modified by applying the mask!!
    data = deepcopy(model_data)
    # only take values where none (!) of the models has infinite values!! (Not just the two that are compared to one another)
    nbModels = length(dims(data, :member))
    maskMissing = dropdims(any(ismissing, data, dims = :member), dims = :member)

    matrixS = zeros(nbModels, nbModels)
    for (i, model_i) in enumerate(eachslice(data[:, :, 1:(end-1)]; dims = :member))
        model_i_mat = Array(model_i)
        model_i_mat[maskMissing .== 1] .= 0
        for (j, model_j) in enumerate(eachslice(data[:, :, (i+1):end]; dims = :member))
            model_j_mat = Array(model_j)
            idx = j + i
            model_j_mat[maskMissing .== 1] .= 0
            s_ij = areaWeightedRMSE(model_i, model_j, maskMissing)
            matrixS[i, idx] = s_ij
        end
    end
    symDistMatrix = matrixS .+ matrixS'
    dim = Array(dims(model_data, :member))
    return YAXArray(
        (Dim{:member1}(dim), Dim{:member2}(dim)), symDistMatrix, deepcopy(model_data.properties),
    )
end


"""
    generalizedDistances(distances_all::YAXArray, weights::DimArray)

For every variable in `distances_all`, compute the weighted sum of all diagnostics.
"""
function generalizedDistances(distances_all::YAXArray, weights::YAXArray)
    model2data = hasdim(distances_all, :member)
    dimensions = model2data ? (:member,) : (:member1, :member2)
    norm = dropdims(median(distances_all, dims=dimensions), dims=dimensions)
    normalized_distances = @d distances_all ./ norm

    distances = model2data ? 
        summarizeEnsembleMembersVector(normalized_distances, false) :
        averageEnsembleMembersMatrix(normalized_distances, false)
    
    weighted_dists = @d distances .* weights
    return dropdims(sum(weighted_dists, dims = :diagnostic), dims=:diagnostic)
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


function getMask(orog_data::YAXArray; mask_out_land::Bool = true)
    ocean_mask = orog_data .== 0
    indices_missing = findall(x -> ismissing(x), ocean_mask)
    # load data into memory with Array() for modification
    ocean_mask_mat = Array(ocean_mask)
    ocean_mask_mat[indices_missing] .= false
    meta = deepcopy(orog_data.properties)
    meta["_ref_id_mask"] = orog_data.properties["_id"]
    mask_arr = YAXArray(dims(ocean_mask), Bool.(ocean_mask_mat), meta)
    return mask_out_land ? mask_arr : mask_arr .== false
end


function getSharedMembers(data::DataMap)
    if !(all(((id, dat),) -> hasdim(dat, :member), data))
        throw(ArgumentError("All datasets must have dimension :member!"))
    end
    members = map(x -> dims(x, :member), collect(values(data)))
    return reduce(intersect, members)
end


function getSharedModels(data::DataMap)
    all_models = Vector(undef, length(data))
    for (i, id) in enumerate(keys(data))
        ds = data[id]
        if hasdim(ds, :model)
            all_models[i] = dims(ds, :model)
        elseif hasdim(ds, :member)
            all_models[i] =
                map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), dims(ds, :member))
        else
            throw(ArgumentError("Data for $id must have dimension :member or :model!"))
        end
    end
    return reduce(intersect, all_models)
end


"""
    filterModels(data::YAXArray, remaining_models::Vector{String})

Remove models from `data` that aren't in `remaining_models` and adapt the 
metadata of `data` accordingly such that 'member_names' and 'model_names' 
corresponds to the remaining data.
"""
function filterModels(data::YAXArray, remaining_models::Vector{String})
    meta = deepcopy(data.properties)
    data =
        hasdim(data, :member) ? data[member=Where(x->x in remaining_models)] :
        data[model=Where(x->x in remaining_models)]

    if haskey(meta, "member_names") && haskey(meta, "model_names")
        indices_out =
            hasdim(data, :member) ?
            findall(x -> !(x in remaining_models), meta["member_names"]) :
            findall(x -> !(x in remaining_models), meta["model_names"])
        deleteat!(meta["model_names"], indices_out)
        deleteat!(meta["member_names"], indices_out)
        data = YAXArray(dims(data), Array(data), meta)
    end
    return data
end


function createGlobalMetaDataDict(
    variable,
    experiment,
    statistic,
    alias,
    timerange;
    paths = Vector{String}(),
)
    metadata = Dict{String,Any}()
    metadata["_paths"] = paths
    metadata["_variable"] = variable
    metadata["_experiment"] = experiment
    metadata["_statistic"] = statistic
    metadata["_alias"] = alias
    metadata["_timerange"] = timerange
    metadata["_id"] = buildMetaDataID(variable, statistic, alias)
    return metadata
end


"""
    loadModelsFromCSV(
        path::String, col_models::String; col_variants::Union{String,Nothing}=nothing
    )

Return a vector with models retrieved from csv file at `path`. If col_variants is provided, 
returned models are on level of model members (MODEL#variant, e.g. AWI-ESM#r1i1p1f1).
"""
function loadModelsFromCSV(
    path::String, col_models::String; col_variants::Union{String, Nothing}=nothing
)
    models = DataFrame(CSV.File(path))
    dropmissing!(models, col_models)
    ids = models[!, col_models]
    if !isnothing(col_variants)
        variants = map(s -> strip.(split(s, ";")), models[!, col_variants])
        for (i, m) in enumerate(models[!, col_models])
            variants[i] = map(v -> join([m, v], MODEL_MEMBER_DELIM), variants[i])
        end
        ids = vcat(variants...)
    end
    return String.(ids)
end


"""
    areaWeightedRMSE(m1::AbstractArray, m2::AbstractArray, mask::AbstractArray)

Compute the area weighted (approximated by cosine of latitudes in radians) root mean squared 
error between `m1` and `m2`. 

# Arguments:
- `m1`: must have dimensions 'lon', 'lat'.
- `m2`: must have dimensions 'lon', 'lat'.
- `mask`: has values 0,1. Locations where mask is 1 get a weight of 0.
"""
function areaWeightedRMSE(m1::AbstractArray, m2::AbstractArray, mask::AbstractArray)
    if dims(m1, :lon) != dims(m2, :lon) || dims(m1, :lat) != dims(m2, :lat)
        msg = "To compute area weigehted RMSE, $m1 and $m2 must be defined on the same lon,lat-grid!"
        throw(ArgumentError(msg))
    end
    squared_diff = (m1 .- m2) .^ 2
    areaweights_mat =
        makeAreaWeightMatrix(Array(dims(m1, :lon)), Array(dims(m1, :lat)); mask)
    weighted_vals = areaweights_mat .* squared_diff
    return sqrt(sum(skipmissing(weighted_vals)))
end


"""
    normalizes(data::Dict{String, Number})

Normalize values for every entry in `data` such that they sum up to 1.

The returned YAXArray has dimension 'diagnostic' whose lookup names are the keys of `data`.
"""
function normalize(data::Dict{String, Number})
    data = filter(((k, v),) -> v != 0, data)
    total = sum(values(data))
    normalized_data = YAXArray(
        (Dim{:diagnostic}(collect(keys(data))),), 
        Array{Float64}(undef, length(data))
    )
    for key in keys(data)
        normalized_data[diagnostic = At(key)] = data[key] / total
    end
    return normalized_data
end



"""
    uncertaintyRanges(
        data::AbstractArray; w::Union{DimArray, Nothing}=nothing, quantiles=[0.167, 0.833]
    )

Compute weighted `quantiles` of timeseries `data`.

# Arguments:
- `data::AbstractArray`: must have dimensions 'time', 'model'.
- `w::Union{DimArray, Nothing}=nothing`: must have dimensions 'model' and 'weight', each 
'weight'-vector must sum up to 1.
- `quantiles=[0.167, 0.833]`: vector with two entries between 0 and 1 representing  the 
lower and upper bound in this order.
"""
function uncertaintyRanges(
    data::AbstractArray; w::Union{YAXArray, Nothing} = nothing, quantiles = [0.167, 0.833]
)
    timesteps = dims(data, :time)
    meta = deepcopy(data.properties)
    meta["_quantiles"] = string.(quantiles)
    uncertainty_ranges = isnothing(w) ?
        YAXArray(
            (dims(data, :time), Dim{:confidence}(["lower", "upper"])),
            Array{Union{Missing, Float64}}(undef, length(timesteps), 2),
            meta
        ) :
        YAXArray(
            (dims(data, :time), dims(w, :weight), Dim{:confidence}(["lower", "upper"])),
            Array{Union{Missing, Float64}}(undef, length(timesteps), length(w.weight), 2),
            meta
        )
    for t in timesteps
        arr = Array(data[time=At(t)])
        if sum(ismissing.(arr)) >= length(arr) - 1 # at least two values must be given
            lower, upper = missing, missing
        else
            values = collect(skipmissing(arr))
            if isnothing(w)
                lower, upper = interpolatedWeightedQuantiles(quantiles, values)
                uncertainty_ranges[time=At(t), confidence=At("lower")] = lower
                uncertainty_ranges[time=At(t), confidence=At("upper")] = upper
            else
                for weight_id in collect(w.weight)
                    lower, upper = interpolatedWeightedQuantiles(
                        quantiles, values; weights = w[weight = At(weight_id)]
                    )
                    uncertainty_ranges[time=At(t), confidence=At("lower"), weight=At(weight_id)] = lower
                    uncertainty_ranges[time=At(t), confidence=At("upper"), weight=At(weight_id)] = upper
                end
            end
        end
    end
    return uncertainty_ranges
end


"""
    interpolatedWeightedQuantiles(
        quantiles::Vector{<:Number}, vals::Vector; weights=nothing
    )

This implementation follows the one used by Brunner et al.
"""
function interpolatedWeightedQuantiles(
    quantiles::Vector{<:Number}, vals::Vector; weights = nothing
)
    if isnothing(weights)
        weights = Array{eltype(vals)}(undef, length(vals))
        weights[ismissing.(vals) .== false] .= 1
    end
    indicesSorted = Array(sortperm(vals)) # gives indices of smallest to largest data point
    weightsSorted = weights[indicesSorted]
    weightedQuantiles = cumsum(weightsSorted) - 0.5 * weightsSorted
    weightedQuantiles = reshape(weightedQuantiles, length(weightedQuantiles), 1)
    # TODO: recheck missing values!
    weightedQuantiles =
        (weightedQuantiles .- minimum(skipmissing(weightedQuantiles))) ./
        maximum(skipmissing(weightedQuantiles))

    interp_linear = Interpolations.linear_interpolation(
        vec(weightedQuantiles),
        vals[indicesSorted],
        extrapolation_bc = Interpolations.Line(),
    )
    return interp_linear(quantiles)
end