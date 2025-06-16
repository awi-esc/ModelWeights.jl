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


function buildCMIP5EnsembleMember(
    realization::Number, initialization::Number, physics::Number,
)
    return "r" * string(realization) * "i" * string(initialization) * "p" * string(physics)
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
    metaAttributesFromESMValToolRecipes(
        base_path_configs::String; constraint::Union{Dict, Nothing} = nothing
    )

Read variable, statistic, experiment and timerange/alias values from ESMValTool recipes 
stored at `base_path_configs` into a vector of Dictionaries storing the respective readoff 
values.

If `constraint` is given, a combination of 'statistic', 'experiment', etc. is only retained 
if each of the respective readoff values aligns with the specified constraints.
"""
function metaAttributesFromESMValToolRecipes(
    base_path_configs::String; constraint::Union{Dict,Nothing} = nothing
)
    paths_configs = filter(
        x -> isfile(x) && endswith(x, ".yml"),
        readdir(base_path_configs, join = true),
    )
    meta_attribs = Vector{MetaData}()
    for path_config in paths_configs
        config = YAML.load_file(path_config)
        data_all = config["diagnostics"]
        aliases = keys(data_all)

        for alias in aliases
            data = data_all[alias]["variables"]
            for (k, v) in data
                #statistic = occursin("_", k) ? String.(split(k, "_"))[2] : ""
                variable = get(v, "short_name", "")
                var_long = get(v, "variable_long_name", "")
                preprocessor = get(v, "preprocessor", "")
                if typeof(v["exp"]) <: String
                    experiment = v["exp"]
                else
                    experiment = join(v["exp"], "-")
                end
                timerange = replace(get(v, "timerange", "full"), "/" => "-")
                meta = createMetaDataFromAttribs(
                    variable, experiment, timerange, k, alias; var_long, preprocessor
                )
                # this may happen for instance if there are different recipes for CMIP5 and CMIP6
                present = any(x -> isEqual(meta, x), meta_attribs)
                if !present
                    push!(meta_attribs, meta)
                end
            end
        end
    end
    if !isnothing(constraint)
        constrainMetaData!(meta_attribs, constraint)
    end
    return meta_attribs
end


function isEqual(m1::MetaData, m2::MetaData)
    return all(x -> getfield(m1, x) == getfield(m2, x), fieldnames(MetaData))
end

function metadataToDict(meta::MetaData; exclude::Union{Vector{Symbol}, Nothing}=nothing)
    d = Dict{String, Any}()
    for f in filter(x -> !(x in exclude), fieldnames(MetaData))
        d[String(f)] = getfield(meta, f)
    end
    return d
end


function renameDictKeys!(data::Dict, keys::Vector)
    for (old_k, new_k) in keys 
        data[new_k] = data[old_k]
        delete!(data, old_k)
    end
end


function alignTimeseries!(data::Vector{YAXArray})
    if !all(map(x -> hasdim(x, :time), data))
        throw(ArgumentError("All datasets must have time dimension to align timeseries!"))
    end
    year_min = minimum(map(x -> minimum(map(Dates.year, dims(x, :time))), data))
    year_max = maximum(map(x -> maximum(map(Dates.year, dims(x, :time))), data))
    nb_years = year_max - year_min + 1

    timerange = DateTime(year_min):Year(1):DateTime(year_max)
    for (i, ds) in enumerate(data)
        s = map(length, otherdims(ds, :time))
            # if ds allows missing values, undef is initialized with missing
        dat = Array{eltype(ds)}(undef, s..., nb_years)
        ds_extended = YAXArray((otherdims(ds, :time)..., Dim{:time}(timerange)), dat)
        ds_extended[time = Where(x->Dates.year(x) in map(Dates.year, dims(ds, :time)))] = ds # ds[time=:]
        data[i] = ds_extended
    end
end


function addPaths!(
    attributes::Vector{MetaData}, 
    path_data::String, 
    dir_per_var::Bool;
    constraint::Union{Dict, Nothing} = nothing 
)
    for meta in attributes
        meta.paths = resolvePaths(meta, path_data, dir_per_var; constraint)
    end
    return nothing
end


# for configuration with yaml file
function get_required_fields_config(ds::Dict)
    data = Dict(
        "base_dir" => get(ds, "base_dir", nothing),
        "exp" => get(ds, "exp", nothing),
        "variables" => get(ds, "variables", nothing),
    )
    if any(isnothing.(values(data)))
        msg = "Config yaml file must specify values for the following keys: $(keys(data))."
        throw(ArgumentError(msg))
    end
    # for fixed variables, no statistics are computed!
    # may add more than 'orog' for which no warning is thrown.
    stats = get(ds, "statistics", "")
    if isempty(stats)
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
    metaDataFromYAML(
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
function metaDataFromYAML(
    content::Dict, is_model_data::Bool; arg_constraint::Union{Dict, Nothing} = nothing
)
    datasets = content["datasets"]
    base_path = get(content, "path_data", "")
    timerange_to_alias = get(content, "timerange_to_alias", Dict{String,String}())

    meta_data = Dict{String, MetaData}()
    for ds in datasets
        meta_ds = Dict{String, MetaData}()
        # get data from config file
        req_fields = get_required_fields_config(ds)
        optional_fields = get_optional_fields_config(ds, timerange_to_alias)
        ds_constraint = merge(req_fields, optional_fields)
        considered_keys = keys(ds_constraint)
        # potentially update with constraint from argument (has precedence over value in yaml file)
        if !isnothing(arg_constraint)
            #setConstraintVal!(ds_constraint, arg_constraint)
            for (field, val) in arg_constraint
                # if subset_shared is given as argument, it will subset considering ALL 
                # loaded datasets, not individual ones!
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
                    meta = createMetaDataFromAttribs(
                        clim_var,
                        ds_constraint["exp"],
                        timerange,
                        clim_var * "_" * stats, # TODO should also be possible to accept any name
                        alias
                    )
                    # for observational data, experiment doesn't make sense
                    if !is_model_data
                        meta["experiment"] = ""
                    end
                    meta.paths = resolvePaths(
                        meta, path_data, ds_constraint["dir_per_var"]; 
                        constraint = ds_constraint
                    )
                    if !isempty(meta.paths)
                        meta_ds[meta.id] = meta
                    else
                        @warn "No data found for $(meta.id)"
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
                meta_data[id].paths = mergeMetaDataPaths(meta_data[id], meta)
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
    constrainFilenames(filenames::Vector{String}, constraint::Dict)

Return subvector of `filenames` retaining those strings that contain at least one of the 
elements in `constraints["models"]` and in `constraints["projects"]`, respectively.
"""
function constrainFilenames(filenames::Vector{String}, constraint::Dict)
    model_constraints = get(constraint, "models", Vector{String}())
    project_constraints = get(constraint, "projects", Vector{String}())
    
    function doIncludeFile(file)
        keep = !isempty(project_constraints) ?
            any([occursin(name, file) for name in project_constraints]) : true
        if keep
            keep = !isempty(model_constraints) ? 
                applyModelConstraints(file, model_constraints) : true
        end
        return keep
    end
    mask = map(f -> doIncludeFile(f), filenames)
    return filenames[mask]
end


function constrainVarSubdirs!(base_paths::Vector{String}, constraint::Dict)
    subdir_constraints = get(constraint, "subdirs", Vector{String}())
    if !isnothing(subdir_constraints) && !isempty(subdir_constraints)
        filter!(p -> any([occursin(name, p) for name in subdir_constraints]), base_paths)
    end
    return nothing
end


"""
    resolvePaths(
        meta::MetaData,
        base_path::String,
        dir_per_var::Bool;
        constraint::Union{Dict, Nothing} = nothing
    )

Return paths to data files for data specified in `meta`, possibly constraint by values in 
`constraint`. The paths were the data is stored is expected to follow the following 
structure (corresponding to the output from ESMValTool used for preprocessing the data):

`base_path` is the top-level directory. If `dir_per_var` is true, `base_path` is assumed to 
have a (or several) subdirectory for each climate variable with _VAR as part of the 
subdirectory's name (e.g. _tas, cmip5_tas, etc.). These subdirectories may be constraint by 
containing at least one of the values in `constraint["subdirs"]`. 

Let BASE refer to `base_path`, or respectively, to the subdirectories for the climate 
variables. Then the following structure is: BASE/preproc/meta.alias/meta.subdir. 
In ESMValTool, `meta.alias` corresponds to the (self-chosen) name under the section 
'diagnostics' and `meta.subdir` to the (self-chosen) name under the section 'variables'.

The returned paths are the paths to all files within this directory, possibly constraint by 
the filenames containig at least one string in `constraint["projects"]` and respectively 
at least one string in `constraint["models"]`. 
"""
function resolvePaths(
    meta::MetaData, 
    base_path::String, 
    dir_per_var::Bool;
    constraint::Union{Dict, Nothing} = nothing
)
    has_constraint = !isnothing(constraint) && !isempty(constraint)
    if dir_per_var
        base_paths = filter(isdir, readdir(base_path, join = true))
        filter!(x -> occursin("_" * meta.variable, x), base_paths)
        if has_constraint
            constrainVarSubdirs!(base_paths, constraint)
        end
    else
        base_paths = [base_path]
    end
    # NOTE: particular data structure assumed here!
    paths_data_dirs = map(base_paths) do p 
        path_data = joinpath(p, "preproc", meta.alias, meta.subdir)
        isdir(path_data) ? path_data : ""
    end
    filter!(x -> !isempty(x), paths_data_dirs)
    if isempty(paths_data_dirs)
        throw(ArgumentError("No directories found at $(base_paths)"))
    end

    paths_to_files = Vector{String}()
    for target_dir in paths_data_dirs
        paths = filter(x -> isfile(x) && endswith(x, ".nc"), readdir(target_dir; join=true))
        paths = has_constraint ? constrainFilenames(paths, constraint) : paths
        append!(paths_to_files, paths)
    end
    return paths_to_files
end


"""
    constrainMetaData!(meta_attributes::Vector{MetaData}, constraint::Dict)

Subset entries in `meta_attributes` so that only those with properties specified 
in `constraint` remain.

# Arguments
- `constraint::Dict`: Mapping to vector specifiying the properties of which at least one 
must be present for an id to be retained.
"""
function constrainMetaData!(meta_attributes::Vector{MetaData}, constraint::Dict)
    timerange_constraints = get(constraint, "timeranges", Vector{String}())
    alias_constraints = get(constraint, "aliases", Vector{String}())
    timerangeOk(attrib::MetaData) = any(x -> attrib.timerange == x, timerange_constraints)
    aliasOk(attrib::MetaData) = any(x -> attrib.alias == x, alias_constraints)

    # timerange and alias don't have to match, it's sufficient if either timerange or alias match
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
    stats_constraints = get(constraint, "statistics", Vector{String}())
    vars_constraints = get(constraint, "variables", Vector{String}())
    stats_ok(attrib::MetaData) = isempty(stats_constraints) || any(x -> attrib.statistic == x, stats_constraints)
    vars_ok(attrib::MetaData) = isempty(vars_constraints) || any(x -> attrib.variable == x, vars_constraints)
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



function distancesData(data::ClimateData, config::Dict{String, Number})
    return distancesData(data.models, data.obs, config)
end

function distancesData(data::ClimateData, diagnostics::Vector{String})
    return distancesData(data.models, data.obs, diagnostics)
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
    # TODO: recheck the metadata here! standard_name should be converted to a vector?!
    distances_all = map(diagnostics) do key 
        distancesData(model_data[key], obs_data[key])
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
    # I think this is not necessary actually (the deepcopy) åWrite test to be sure!!
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
    distances_all = map(diagnostics) do key 
        distancesModels(model_data[key])
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


function fixModelNameMetadata(name::String)
    return get(MODEL_NAME_FIXES, name, name)
end


function getMask(orog_data::YAXArray; mask_out_land::Bool = true)
    ocean_mask = orog_data .== 0
    indices_missing = findall(x -> ismissing(x), ocean_mask)
    # load data into memory with Array() for modification
    ocean_mask_mat = Array(ocean_mask)
    ocean_mask_mat[indices_missing] .= false
    meta = deepcopy(orog_data.properties)
    meta["_ref_id_mask"] = orog_data.properties["id"]
    mask_arr = YAXArray(dims(ocean_mask), Bool.(ocean_mask_mat), meta)
    return mask_out_land ? mask_arr : mask_arr .== false
end


function sharedMembers(data::DataMap)
    if !(all(((id, dat),) -> hasdim(dat, :member), data))
        throw(ArgumentError("All datasets must have dimension :member!"))
    end
    members = map(x -> dims(x, :member), collect(values(data)))
    return reduce(intersect, members)
end


function sharedModels(data::DataMap)
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


function metaVecToDict(attributes::Vector{MetaData})
    datasets_meta = Dict{String, MetaData}()
    for meta in attributes
        if !isempty(meta.paths)
            id = meta.id
            if haskey(datasets_meta, id)
                @warn "several times same data id for $id"
                datasets_meta[id].paths = mergeMetaDataPaths(datasets_meta[id], meta)
            else
                datasets_meta[id] = meta
            end
        end
    end
    return datasets_meta
end


function createMetaDataFromAttribs(
    var_short::String, 
    experiment::String, 
    timerange::String, 
    subdir::String, 
    alias::String; 
    var_long::String = "", 
    preprocessor::String = "", 
    paths = Vector{String}()
)
    var_stat = String.(split(subdir, "_"))
    return MetaData(
        paths = paths, 
        variable = var_short, 
        experiment = experiment,
        timerange = timerange,
        alias = alias,
        subdir = subdir,
        variable_long = var_long, 
        statistic = length(var_stat) == 2 ? var_stat[2] : "", # TODO
        esmvaltoolPreproc = preprocessor,
        id = buildMetaDataID(subdir, alias)
    )
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
    normalize(data::Dict{String, <:Number})

Normalize values for every entry in `data` such that they sum up to 1. If remove_zero is
true (default), the returned dictionary does not contain entries for which values were 0.
"""
function normalize(data::Dict{String, <:Number}; remove_zero::Bool=true)
    result = Dict{String, Float64}()
    total = sum(values(data))
    data = remove_zero ? filter(((k, v),) -> v != 0, data) : data
    for k in keys(data)
        result[k] = data[k] ./ total 
    end
    return result
end


"""
    normalizeToYAX(data::Dict{String, <:Number})

Normalize values for every entry in `data` such that they sum up to 1 and return a YAXArray 
with dimension 'diagnostic' whose lookup names are the keys of `data`.
"""
function normalizeToYAX(data::Dict{String, <:Number})
    normed = normalize(data)
    normalized_yax = YAXArray(
        (Dim{:diagnostic}(collect(keys(normed))),), Array{Float64}(undef, length(normed))
    )
    for key in keys(normed)
        normalized_yax[diagnostic = At(key)] = normed[key]
    end
    return normalized_yax
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