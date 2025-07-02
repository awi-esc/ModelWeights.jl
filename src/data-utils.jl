"""
    buildMemberID(meta::Dict, model::String)

Return member id for `model` which is the model name itself followed by '#' and the variant 
label (simulation id).

The unique member ids correspond to the variant labels of CMIP6 models, e.g. r1i1p1f1.
For CMIP5 models this id is built from the respective metadata.
For CMIP6 models, the grid is further added to the end of the id seperated with _.

# Arguments:
- `meta::Dict`: For CMIP5 data, must have keys: 'mip_era', 'realization', 
'initialization_method', 'physics_version'. For CMIP6 data must have keys: 'variant_label', 
'grid_label'.
"""
function buildMemberID(meta::Dict, model::String)
    fn_err(k::String) = throw(ArgumentError("$k not defined in metadata!"))
    mip_era = get(() -> fn_err("mip_era"), meta, "mip_era")
    member = ""
    if lowercase(mip_era) == "cmip5"
        variant = buildCMIP5EnsembleMember(
            get(() -> fn_err("realization"), meta, "realization"),
            get(() -> fn_err("initialization_method"), meta, "initialization_method"),
            get(() -> fn_err("physics_version"), meta, "physics_version")
        )
        member = join([model, variant], MODEL_MEMBER_DELIM)
    elseif lowercase(mip_era) == "cmip6"
        variant = get(() -> fn_err("variant_label"), meta, "variant_label")
        grid = get(() -> fn_err("grid_label"), meta, "grid_label")
        member = join([model, variant, grid], MODEL_MEMBER_DELIM, "_")
    else 
        throw(ArgumentError("Only CMIP5 and CMIP6 supported so far."))
    end
    return member
end


"""
    subsetModelData(data::YAXArray, shared_models::Vector{String})

Return subset of `data` containing only data from models in `shared_models`. 

Takes care of metadata.

# Arguments:
- `data`: must have dimension 'member' or 'model'
- `shared_models`: models, which can either be on level of models or members of models 
('modelname#memberID[_grid]').
"""
function subsetModelData(data::YAXArray, shared_models::Vector{String})
    if isempty(shared_models)
        @warn "Vector of models to subset data to is empty!"
        return data
    end
    data = deepcopy(data)
    dim_symbol = hasdim(data, :member) ? :member : :model
    dim_names = Array(dims(data, dim_symbol))
    if dim_symbol == :member
        models = map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), shared_models)
        # if shared_models is on the level of models, the following should be empty
        # otherwise, nothing is filtered out, and members is the same as shared_models 
        members = filter(x -> !(x in models), shared_models)
        if !isempty(members) # shared models on level of members
            indices = findall(m -> m in members, dim_names)
        else
            # important not to use dim_names here, since e.g. model=AWI would be found in dim_names where model is actually AWI-X for instance
            models_data = modelsFromMemberIDs(dim_names) # NOTE: should yield same: models_data = data.properties["model_names"]
            indices = findall(m -> m in models, models_data)
        end
        data = data[member=indices]
    else
        indices = findall(m -> m in shared_models, dim_names)
        data = data[model=indices]
    end
    # also subset Metadata vectors!
    attributes = filter(k -> data.properties[k] isa Vector, keys(data.properties))
    for key in attributes
        data.properties[key] = data.properties[key][indices]
    end
    return data
end


"""
    mergeMetaDataFromMultipleFiles(data::Vector{<:YAXArray})

Combine arrays in `data` into a single YAXArray with meta combined from all datasets into 
lists, with missing if key wasnt in a dataset.
"""
function mergeMetaDataFromMultipleFiles(data::Vector{<:YAXArray})
    n_files = length(data)
    meta_keys = unique(vcat(map(x -> collect(keys(x.properties)), data)...))
    meta_dict = Dict{String, Any}()
    for (i, ds) in enumerate(data)
        for key in meta_keys
            values = get!(meta_dict, key, repeat(Any[missing], outer=n_files))
            values[i] = get(ds.properties, key, missing)
        end
    end
    return meta_dict
end


"""
    combineModelsFromMultipleFiles(
        data::Vector{<:YAXArray}; 
        names::Vector{String} = Vector{String}(),
        meta::Union{Dict{String, T}, Nothing} = nothing,
        dtype::String = "undef",
        new_dim::Symbol = :model,
        sorted::Bool = true
) where T <: Any

Combine `data` from different files into a single YAXArray. 
The meta data of the returned YAXArray is a combination of all values from all datasets 
using vectors with missing values for datasets that didn't have the respective metadata 
entry in their metadata, further combined with the properties in `meta` if provided.

All elements in `data` must share all dimensions except for time (if present).
For timeseries data, the time dimension may cover different ranges. In that case the 
maximal timeseries is used and filled with NaN for missing values.  

The combined YAXArray has the additional dimension `new_dim` (default: :model) with `names`
as values. If `names` is not provided, default values 'model1', 'model2, etc. are used 
instead.

If `sorted` is true, model dimension of returned data is sorted alphabetically and the 
vector entries in the metadata dictionary of the returned array are also sorted accordingly.
"""
function combineModelsFromMultipleFiles(
    data::Vector{<:YAXArray}; 
    names::Vector{String} = Vector{String}(),
    meta::Union{Dict{String, T}, Nothing} = nothing,
    new_dim::Symbol = :model,
    sorted::Bool = true
) where T <: Any
    if isempty(data)
        @warn "Data vector is empty!"
        return nothing
    end
    data_sizes = unique(map(size, data))
    if length(data_sizes) != 1
        if !all(map(x -> hasdim(x, :time), data))
            msg = "Data does not have the same size across all models: $(data_sizes)"
            throw(ArgumentError(msg))
        else
            # if difference only in time, use maximal possible timeseries and add NaNs
            alignTimeseries!(data)
        end
    end
    use_default_names = isempty(names)
    if use_default_names
        names = map(x -> string(new_dim) * string(x), 1:length(data))
    end
    if sorted && !use_default_names
        sort_indices = sortperm(names)
        names = names[sort_indices]
        data = data[sort_indices]
    end
    meta_dict = mergeMetaDataFromMultipleFiles(data)    
    meta_dict["info"] = isnothing(meta) ? Dict() : meta
    dimData = concatenatecubes(data, Dim{new_dim}(names))
    dimData = YAXArray(dimData.axes, dimData.data, meta_dict)
    if length(names) != length(unique(names))
        duplicates = unique([m for m in names if sum(names .== m) > 1])
        @warn "Some datasets appear more than once" duplicates
    end
    return dimData
end


"""
    memberIDsFromPaths(all_paths::Vector{String})

For every path in `all_paths` return a string of the form modelname#memberID[_grid]
that identifies the corresponding model member. Filenames must follow the CMIP-standard
(proj_name_mip_exp_id_variable[_grid].nc).

For CMIP6 models the abbreviation of the grid is added to the model name.
"""
function memberIDsFromPaths(all_paths::Vector{String})
    all_filenames = split.(basename.(all_paths), "_")
    all_members = Vector{String}(undef, length(all_filenames))
    for (i, fn_parts) in enumerate(all_filenames)
        model = join(fn_parts[[2, 5]], MODEL_MEMBER_DELIM)
        # add grid to model name for CMIP6 models:
        if fn_parts[1] != "CMIP5"
            model = splitext(model * "_" * fn_parts[7])[1]
        end
        all_members[i] = model
    end
    return unique(all_members)
end


"""
    memberIDFromFilenameMeta(fn_meta::FilenameMeta)

Return a string of the form modelname#memberID[_grid] that identifies the corresponding model member.

For CMIP6 models the abbreviation of the grid is added to the model name.
"""
function memberIDFromFilenameMeta(fn_meta::FilenameMeta, mip_era::String)
    id = join([fn_meta.model, fn_meta.variant], MODEL_MEMBER_DELIM)
    if lowercase(mip_era) == "cmip6"
        id = join([id, fn_meta.grid], "_")
    end
    return id
end


"""
    findModelInPaths(model_id::String, paths::Vector{String})

Return true if `model_id` is found in filename of any path in `paths`, else false.

The paths are assumed to follow the standard CMIP filename structure, i.e. <variable_id>_<table_id>_<source_id>_<experiment_id >_<member_id>_<grid_label>[_<time_range>].nc([see here](https://docs.google.com/document/d/1h0r8RZr_f3-8egBMMh7aqLwy3snpD6_MrDz1q8n5XUk/edit?tab=t.0)).

# Arguments:
- `model_id::String`: has form modelname[#memberID[_grid]]
- `paths::Vector{String}`: paths to be searched
"""
function findModelInPaths(model_id::String, paths::Vector{String})
    model_parts = String.(split(model_id, MODEL_MEMBER_DELIM))
    model = model_parts[1]
    has_member = length(model_parts) == 2
    member_grid = has_member ? split(model_parts[2], "_") : nothing
    has_grid = !isnothing(member_grid) && length(member_grid) == 2
    member = has_member ? member_grid[1] : nothing
    grid = has_grid ? member_grid[2] : nothing

    is_found = false
    filenames = map(basename, paths)
    for fn in filenames
        found_model = occursin("_" * model * "_", fn)
        if !found_model
            continue
        else
            if has_member
                found_member = occursin("_" * member * "_", fn)
                if found_member
                    if has_grid
                        fn_no_ending = splitext(fn)[1]
                        found_grid = endswith(fn_no_ending, "_" * grid) || occursin("_" * grid * "_", fn)
                        if found_grid
                            is_found = true
                            break
                        else
                            continue
                        end
                    else
                        is_found = true
                        break
                    end
                else # member not found -> continue with next path
                    continue
                end
            else
                is_found = true
                break
            end
        end
    end
    return is_found
end


"""
    getCMIPModelsKey(meta::Dict)

Return the respective key to retrieve model names in CMIP6 ('source_id') and 
CMIP5 ('model_id') data.

If both keys are present, 'source_id' used in CMIP6 models is returned, if none 
is present, throw ArgumentError.
"""
function getCMIPModelsKey(meta::Dict)
    attributes = keys(meta)
    if "source_id" in attributes
        if "model_id" in attributes
            msg1 = "Dictionary contains keys source_id and model_id, source_id is used! "
            @debug msg1 meta["source_id"] meta["model_id"]
        end
        return "source_id"
    elseif "model_id" in attributes
        return "model_id"
    else
        msg = "Metadata must contain one of 'source_id' (pointing to names of CMIP6 models) or 'model_id' (CMIP5)."
        throw(ArgumentError(msg))
    end
end


"""
    filterPathsSharedModels(paths::Vector{String}, shared_models::Vector{String})

Every vector of paths in `all_paths` is filtered s.t. it only contains models or model 
members given in `shared_models`.
"""
function filterPathsSharedModels(paths::Vector{String}, shared_models::Vector{String})
    if isempty(shared_models)
        @warn "No models shared across data!"
        return nothing
    end
    mask = map(p -> applyModelConstraints(p, shared_models), paths)
    return paths[mask]
end

function filterPathsSharedModels(all_paths::Vector{Vector{String}}, level_shared::Level)
    shared = sharedModels(all_paths, level_shared)
    return map(paths -> filterPathsSharedModels(paths, shared), all_paths)
end

function filterPathsSharedModels(
    all_paths::Vector{Vector{String}}, level_shared::Union{String, Symbol}
)
    return filterPathsSharedModels(all_paths, getLevel(level_shared))
end


function modelsFromMemberIDs(members::AbstractVector{<:String}; uniq::Bool=false)
    models = map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), members)
    return uniq ? unique(models) : models
end


"""
    approxAreaWeights(latitudes::Vector{<:Number})

Create a YAXArray with the cosine of `latitudes` which approximates the cell 
area on the respective latitude.
"""
function approxAreaWeights(latitudes::Vector{<:Number})
    # cosine of the latitudes as proxy for grid cell area
    area_weights = cos.(deg2rad.(latitudes))
    return YAXArray((Dim{:lat}(latitudes),), area_weights)
end


function makeAreaWeightMatrix(
    longitudes::Vector{<:Number},
    latitudes::Vector{<:Number};
    mask::Union{AbstractArray,Nothing} = nothing,
)
    area_weights = approxAreaWeights(latitudes)
    area_weighted_mat = repeat(area_weights', length(longitudes), 1)
    if !isnothing(mask)
        area_weighted_mat = ifelse.(mask .== 1, 0, area_weighted_mat)
    end
    area_weighted_mat = area_weighted_mat ./ sum(area_weighted_mat)
    return YAXArray((Dim{:lon}(longitudes), Dim{:lat}(latitudes)), area_weighted_mat)
end

function checkDataStructure(path_data::String, dir_per_var::Bool)
    if !dir_per_var
        subdirs = filter(x -> isdir(joinpath(path_data, x)), readdir(path_data))
        if !("preproc" in subdirs)
            throw(ArgumentError("If dir_per_var is false, path_data must contain a directory named 'preproc'"))
        end
    end
    return nothing
end


function buildMetaDataID(meta::Dict)
    fn_err(k::String) = throw(ArgumentError("$k not defined in metadata!"))
    subdir = get(meta, "subdir", nothing)
    if isnothing(subdir)
        subdir = join(
            [get(() -> fn_err("subdir, variable"), meta, "variable"),
             get(() -> fn_err("subdir, statistic"), meta, "statistic")], 
             "_"
        )
    end
    return join([subdir, get(() -> fn_err("alias"), meta, "alias")], "_")
end


function buildMetaDataID(meta::MetaData)
    subdir = meta.subdir
    if isempty(subdir)
        subdir = (!isempty(meta.statistic) && !isempty(meta.variable)) ? 
            join([meta.variable, meta.statistic], "_") :
            throw(ArgumentError("Neither subdir nor statistic and variable are in metadata!"))
    end
    if isempty(meta.alias)
        throw(ArgumentError("alias (name of diagnostic in ESMValTool) cannot be empty in MetaData!"))
    end
    return join([subdir, meta.alias], "_")
end


"""
    joinDicts(v::Dict...; warn_msg::String="")
"""
function joinDicts(v::Dict...; warn_msg::String="")
    result = typeof(v[1])()
    for dm in v
        shared_keys = sharedKeys(result, dm)
        if length(shared_keys) > 0
            dups = filter(!isempty, map(shared_keys) do x
                    dm[x] != result[x] ? (x, result[x], dm[x]) : ()
                end
            )
            if !isempty(dups) && !isempty(warn_msg)
                @warn warn_msg dups
            end
        end
        result = merge(result, dm)
    end
    return result
end


function physicsFromMember(member::String)
    regex = r"(p\d+)(f\d+)?(_.*)?$"
    return String.(match(regex, member).captures[1])
end


"""
    averageEnsembleMembersMatrix(data::AbstractArray, updateMeta::Bool)

Compute the average across all members of each model for each given variable 
for model to model data, e.g. distances between model pairs.

# Arguments:
- `data`: with at least dimensions 'member1', 'member2'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
Set to false if vectors refer to different variables for instance. 
"""
function averageEnsembleMembersMatrix(data::YAXArray, updateMeta::Bool)
    data = setLookupsFromMemberToModel(data, ["member1", "member2"])
    models = String.(collect(unique(dims(data, :model1))))

    grouped = groupby(data, :model2 => identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims = (:model2,)), grouped)
    combined = cat(averages..., dims = (Dim{:model2}(models)))

    grouped = groupby(combined, :model1 => identity)
    averages = map(entry -> mapslices(Statistics.mean, entry, dims = (:model1,)), grouped)
    combined = cat(averages..., dims = (Dim{:model1}(models)))

    for m in models
        combined[model1=At(m), model2=At(m)] .= 0
    end

    meta = updateMeta ? updateGroupedDataMetadata(data.properties, grouped) : data.properties
    combined = rebuild(combined; metadata = meta)

    l = Lookups.Categorical(sort(models); order = Lookups.ForwardOrdered())
    combined = combined[model1=At(sort(models)), model2=At(sort(models))]
    combined = DimensionalData.Lookups.set(combined, model1 = l, model2 = l)
    return combined
end


function joinDataMaps(v::Union{Dict{String, MetaData}, DataMap}...; warn_msg::String="")
    return joinDicts(v...; warn_msg)
end


"""
    writeDataToDisk(data, target_path::String)

Save `data` as Julia obj if `target_path` has ending '.jld2', otherwise save as binary.
If file at `target_path` already exists, timestamp is added if `overwrite` is false (default).
"""
function writeDataToDisk(data, target_path::String; overwrite::Bool = false)
    target_path = overwrite ? target_path : Data.individuatePath(target_path)
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
    buildCMIP5EnsembleMember(realization::Number, initialization::number, physics::Number)

Build rip-abbreviations from `realization`, `initialization` and `physics`; especially for 
CMIP5 models which do not have it in their metadata.
"""
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
each vector only contains N entries (N=number of models (i.e. without unique members)).

If the metadata for members of a model differ across members, the respective
entry in the vector will be a vector itself.
"""
function updateGroupedDataMetadata(
    meta::Dict, grouped_data::DimensionalData.DimGroupByArray,
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
    metaDataFromESMValToolRecipes(
        base_path_configs::String; constraint::Union{Dict, Nothing} = nothing
    )

Read variable, statistic, experiment and timerange/alias values from ESMValTool recipes 
stored at `base_path_configs` into a vector of Dictionaries storing the respective readoff 
values.
"""
function metaDataFromESMValToolRecipes(
    base_path_configs::String;
    constraint::Union{Dict, Nothing} = nothing 
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
                variable = get(v, "short_name", "")
                var_long = get(v, "variable_long_name", "")
                preprocessor = get(v, "preprocessor", "")
                
                var_stat = String.(split(k, "_"))
                statistic = length(var_stat) == 2 && var_stat[1] == variable ? var_stat[2] : nothing
                if typeof(v["exp"]) <: String
                    experiment = v["exp"]
                else
                    # TODO: when would this happen??
                    experiment = join(v["exp"], "-")
                end
                timerange = replace(get(v, "timerange", "full"), "/" => "-")
                meta = MetaData(
                    variable, experiment, alias; 
                    timerange = timerange,
                    subdir = k,
                    variable_long = var_long,
                    statistic = statistic,
                    esmvaltoolPreproc = preprocessor
                )
                # if there are different recipes for CMIP5 and CMIP6, 'meta' might already 
                # have been added to 'meta_attribs'
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


function metadataToDict(meta::MetaData; exclude::Vector{Symbol}=Vector{Symbol}())
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


function alignTimeseries!(data::Vector{<:YAXArray})
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
        ds_extended = YAXArray((otherdims(ds, :time)..., Dim{:time}(timerange)), dat, ds.properties)
        ds_extended[time = Where(x->Dates.year(x) in map(Dates.year, dims(ds, :time)))] = ds # ds[time=:]
        data[i] = ds_extended
    end
end


"""
    metaDataFromYAML(content::Dict)

Return metadata of data specified in `content` possibly constrained by values in `arg_constraint`.

For constraints that are specified in `content` as well as in the `arg_constraint` argument, 
the values of the latter have precedence over the former. The constraints given
in the argument `arg_constraint` are applied to EVERY dataset specified in the config 
file.

# Arguments:
- `content`: content of config yaml file specifying meta attributes and paths of data
"""
function metaDataFromYAML(ds::Dict)
    fn_err(x) = throw(ArgumentError("$(x) must be provided for each dataset in config yaml file!"))
    experiment = get(() -> fn_err("exp"), ds, "exp")
    variables = get(() -> fn_err("variables"), ds, "variables")
    aliases = get(() -> fn_err("aliases"), ds, "aliases")
    
    subdirs = get(ds, "subdirs", nothing)
    statistics = get(ds, "statistics", nothing)

    meta_ds = Vector{MetaData}()
    if isnothing(subdirs)
        subdirs = isnothing(statistics) ? 
            fn_err("if subdirs not provided, statistics") : 
            combineAll(variables, statistics)
        for alias in aliases
            for subdir in subdirs
                var, stats = string.(split(subdir, "_"))
                push!(meta_ds, MetaData(var, experiment, alias; subdir=subdir, statistic=stats))     
            end
        end
    else
        for var in variables
            for alias in aliases
                for subdir in subdirs
                    push!(meta_ds, MetaData(var, experimen, alias; subdir=subdir))
                end
            end
        end
    end
    return meta_ds
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


function constrainSubdirs!(paths::Vector{String}, subdir_constraints::Vector{String})
    if !isempty(subdir_constraints)
        filter!(p -> any([occursin(name, p) for name in subdir_constraints]), paths)
    end
    return nothing
end


"""
    resolvePathsFromMetaData( 
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
containing at least one of the values in `constraint["base_subdirs"]`. 

Let BASE refer to `base_path`, or respectively, to the subdirectories for the climate 
variables. Then the following structure is: BASE/preproc/meta.alias/meta.subdir. 
In ESMValTool, `meta.alias` corresponds to the (self-chosen) name under the section 
'diagnostics' and `meta.subdir` to the (self-chosen) name under the section 'variables'.

The returned paths are the paths to all files within this directory, possibly constraint by 
the filenames containig at least one string in `constraint["projects"]` and respectively 
at least one string in `constraint["models"]`. 
"""
function resolvePathsFromMetaData(
    meta::MetaData, 
    base_path::String, 
    dir_per_var::Bool;
    constraint::Union{Dict, Nothing} = nothing
)
    has_constraint = !isnothing(constraint) && !isempty(constraint)
    if dir_per_var
        base_paths = filter(isdir, readdir(base_path, join = true)) # all subdirectories
        filter!(x -> occursin("_" * meta.variable, x), base_paths) # just subdirs for variable
        if has_constraint
            constrainSubdirs!(base_paths, get(constraint, "base_subdirs", Vector{String}()))
        end
    else
        base_paths = [base_path]
    end
    paths_data_dirs = map(base_paths) do p 
        path_data = joinpath(p, "preproc", meta.alias, meta.subdir)
        isdir(path_data) ? path_data : ""
    end
    filter!(x -> !isempty(x), paths_data_dirs)
    if isempty(paths_data_dirs)
        throw(ArgumentError("No directories found at $(base_paths)"))
    end
    paths_to_files = map(p -> collectNCFilePaths(p; constraint), paths_data_dirs)
    return vcat(paths_to_files...)
end


function parseFilename(filename::String, format::String)
    # TODO: parse individual format
end


function parseFilename(filename::String, format::Symbol)
    err_msg = "Only filename formats $(keys(FILENAME_FORMATS)) defined. Found: $(format)."
    attribs = get(() -> throw(ArgumentError(err_msg)), FILENAME_FORMATS, format)
    values = Vector{Union{String, Nothing}}(split(filename, "_"))
    diff = length(attribs) - length(values)
    if diff > 0
        append!(values, repeat([nothing], diff))
    elseif diff < 0
        throw(ArgumentError("$(filename) has more parts than expected for format $format!"))
    end
    mapping = Dict(attribs .=> values)
    return FilenameMeta(
        variable = mapping["variable"],
        table_id = mapping["table_id"],
        model =  mapping["model"],
        exp = mapping["exp"],
        variant = mapping["variant"],
        fn = filename,
        grid = get(mapping, "grid", nothing),
        timerange = get(mapping, "timerange", nothing),
        mip = get(mapping, "mip", nothing)
    )
end

function collectNCFilePaths(path_data_dir::String; constraint::Union{Dict, Nothing}=nothing)
    paths_to_files = filter(x -> isfile(x) && endswith(x, ".nc"), readdir(path_data_dir; join=true))
    if isempty(paths_to_files)
        @warn "No .nc files found!"
    end
    if !isempty(paths_to_files) && !isnothing(constraint) && !isempty(constraint)
        paths_to_files = constrainFilenames(paths_to_files, constraint)
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
        keep_file = findModelInPaths(model, [path_model_data])
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
    distances_all = Vector{<:YAXArray}(undef, length(diagnostics))
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


"""
    sharedLevelMembers(data::DataMap)

Return vector of model members that are shared across all entries in `data`.
"""
function sharedLevelMembers(data::DataMap)
    if !(all(((id, dat),) -> hasdim(dat, :member), data))
        throw(ArgumentError("All datasets must have dimension :member!"))
    end
    members = map(x -> dims(x, :member), collect(values(data)))
    return reduce(intersect, members)
end


"""
    sharedLevelModels(data::DataMap)

Return vector of models (on level of models, not members) that are shared across all entries in `data`.
"""
function sharedLevelModels(data::DataMap)
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


"""
    metaDataChecksCMIP(meta::Dict, path::String)

Check model names as retrieved from the metadata for potential inconsistencies wrt filename.
"""
function metaDataChecksCMIP(meta::NCDatasets.CommonDataModel.Attributes, path::String)
    name = meta[getCMIPModelsKey(Dict(meta))]
    if !occursin(name, basename(path)) && !(name in keys(MODEL_NAME_FIXES))
        @warn "model name as read from metadata of stored .nc file ($name) and used as dimension name is not identical to name appearing in filename $(basename(path))."
    end 
    return nothing
end