"""
    subsetModelData(datamap::DataMap; level::Level=MEMBER_LEVEL)

For those datasets in `datamap` that specify data on the level `level` (i.e. have dimension 
:member or :model), return a new DataMap with subset of data s.t. the new datasets all have 
the same models (MODEL_LEVEL) or members (MEMBER_LEVEL).

If no models are shared across datasets, return the input `datamap`.
"""
function subsetModelData(datamap::DataMap; level::Level = MEMBER_LEVEL)
    shared_models = sharedModels(datamap, level)
    if isempty(shared_models)
        @warn "no shared models in datamap"
        return datamap
    end
    subset = DataMap()
    for (id, data) in datamap
        subset[id] = subsetModelData(deepcopy(data), shared_models)
    end
    return subset
end


function sharedModels(data::DataMap, level::Level)
    return level == MEMBER_LEVEL ? sharedLevelMembers(data) : sharedLevelModels(data)
end


function sharedModels(all_paths::Vector{Vector{String}}, subset_shared::Level)
    all_models = memberIDsFromPaths(vcat(all_paths...))
    if subset_shared == MODEL_LEVEL
        all_models = unique(modelsFromMemberIDs(all_models))
    end
    # check whether models are found in paths for every dataset
    indices_shared = []
    for (idx, model) in enumerate(all_models)
        is_found = false
        for paths in all_paths
            is_found = findModelInPaths(model, paths)
            if !is_found
                break
            end
        end
        if is_found
            push!(indices_shared, idx)
        end
    end
    return all_models[indices_shared]
end


"""
    alignPhysics(
        datamap::DataMap, members::Vector{String}; subset_shared::Union{Level, Nothing}=nothing,
    )

Return new DataMap with only the models retained that share the same physics as 
the respective model's members in `members`.

If `subset_shared` is set, resulting DataMap is subset accordingly.
"""
function alignPhysics(
    datamap::DataMap, members::Vector{String}; subset_shared::Union{Level, Nothing}=nothing,
)
    data = deepcopy(datamap)
    models = unique(modelsFromMemberIDs(members))
    for model in models
        # retrieve allowed physics as given in members 
        member_ids = filter(m -> startswith(m, model * MODEL_MEMBER_DELIM), members)
        physics = physicsFromMember.(member_ids)
        for (_, ds) in data
            # filter data s.t. of current model only members with retrieved physics are kept
            model_indices = findall(
                x -> startswith(x, model * MODEL_MEMBER_DELIM),
                Array(dims(ds, :member)),
            )
            indices_out =
                filter(x -> !(ds.properties["physics"][x] in physics), model_indices)
            if !isempty(indices_out)
                indices_keep = filter(x -> !(x in indices_out), 1:length(dims(ds, :member)))
                members_kept = ds.properties["member_names"][indices_keep]
                data[ds.properties["id"]] = subsetModelData(ds, members_kept)
            end
        end
    end
    if !isnothing(subset_shared)
        shared_models = sharedModels(data, subset_shared)
        for (id, model_data) in data
            data[id] = subsetModelData(model_data, shared_models)
        end
    end
    return data
end


""" 
    summarizeEnsembleMembersVector(
        data::YAXArray, updateMeta::Bool; fn::Function=Statistics.mean
)

For each model compute a summary statistic (default: mean) across all its members. 
The returned YAXArray has dimension 'model'.

# Arguments:
- `data::YAXArray`: YAXArray with at least dimension 'member'
- `updateMeta::Bool`: set true if the vectors in the metadata refer to 
different models. Set to false if vectors refer to different variables.
- `fn::Function`: Function to be applied on data
"""
function summarizeEnsembleMembersVector(
    data::YAXArray, updateMeta::Bool; fn::Function = Statistics.mean
)
    throwErrorIfModelDimMissing(data)
    data = setLookupsFromMemberToModel(data, ["member"])
    models = unique(Array(dims(data, :model)))
    dimensions = otherdims(data, :model)
    s = isempty(dimensions) ? (length(models),) : (size(dimensions)..., length(models))
    summarized_data = YAXArray(
        (otherdims(data, :model)..., Dim{:model}(models)),
        Array{eltype(data)}(undef, s),
        deepcopy(data.properties),
    )
    for m in models
        dat = data[model = Where(x -> x == m)]
        average =
            isempty(dimensions) ? fn(dat) : mapslices(x -> fn(x), dat; dims = (:model,))
        summarized_data[model = At(m)] = average
    end
    summarized_data = replace(summarized_data, NaN => missing)

    meta = deepcopy(data.properties)
    # TODO: fix updateMetadat
    #meta = updateMeta ? updateGroupedDataMetadata(meta, grouped) : meta
    return YAXArray(dims(summarized_data), summarized_data.data, meta)
end


"""
    summarizeEnsembleMembersVector!(data::DataMap)

Set values for every dataset in `data` to the average across all members of 
each model.
"""
function  summarizeEnsembleMembersVector!(data::DataMap)
    for (k, ds) in data
        if hasdim(ds, :member)
            @info "average ensemble members for $k"
            data[k] = summarizeEnsembleMembersVector(ds, true)
        end
    end
    return nothing
end


function addMasks!(datamap::DataMap, id_orog_data::String)
    orog_data = datamap[id_orog_data]
    datamap["mask_land"] = getMask(orog_data; mask_out_land = false)
    datamap["mask_ocean"] = getMask(orog_data; mask_out_land = true)
    return nothing
end

### ----------------------------------------------------------------------------------------
###                                LOADING DATA                                           
### ----------------------------------------------------------------------------------------
"""
    loadPreprocData(
        paths::Vector{String},
        fn_format::Symbol;
        sorted::Bool = true, 
        dtype::DataType = MODEL_OBS_DATA, 
        meta::Union{Dict{String, T}, Nothing} = nothing
    )

Return data loaded from `paths` as single YAXArray. The names of the data files are assumed 
to follow the CMIP-standard. 

Filenames that do not start with 'cmip5' or 'cmip6' (upper/lower case doesnt matter), are 
considered observational data.
"""
function loadPreprocData(
    paths::Vector{String},
    fn_format::Symbol;
    sorted::Bool = true, 
    dtype::DataType = MODEL_OBS_DATA, 
    meta_info::Union{Dict{String, T}, Nothing} = nothing
) where T <: Any
    data = Vector{YAXArray}()
    filenames = first.(splitext.(basename.(paths)))
    filenames_meta = parseFilename.(filenames, fn_format)
    for (path, filename, meta) in zip(paths, filenames, filenames_meta)
        @debug "processing file.." * path
        ds = NCDataset(path)
        # check if data is model data or observational data
        mip = get(ds.attrib, "project_id", nothing) # for CMIP5
        mip = isnothing(mip) ? get(ds.attrib, "mip_era", nothing) : mip # for CMIP6
        is_model_data = !isnothing(mip)
        
        @debug "is_model_data: $is_model_data, fn_format: $fn_format, meta:" meta
        clim_var = (is_model_data || fn_format == :esmvaltool) ? 
            meta.variable : get(meta_info, "variable", nothing)
        if isnothing(clim_var)
            throw(ArgumentError("for observational data, variable must be provided in meta_info if filename format is not based on standard assumed by esmvaltool."))
        end

        dsVar = nothing
        try 
            dsVar = ds[clim_var]
        catch e
            @error "Variable '$clim_var' according to filename format '$fn_format' not found!"
            throw(e)
        end
        # if clim_var == "amoc"
        #     dsVar = ds["msftmz"]
        # end
        attributes = merge(dsVar.attrib, ds.attrib)
        if clim_var == "msftmz"
            sector = get(ds, "sector", nothing)
            if !isnothing(sector)
                merge!(attributes, sector.attrib)
            end
        end
        props = Dict{String, Any}() # metadata just for this file!
        props["path"] = path

        if is_model_data && dtype in [MODEL_DATA, MODEL_OBS_DATA] 
            model_key = getCMIPModelsKey(Dict(ds.attrib))
            # add mip_era for models since it is not provided in CMIP5-models
            get!(attributes, "mip_era", "CMIP5")
            # check model names as retrieved from the metadata for potential inconsistencies wrt filename
            name = ds.attrib[model_key] # just the model name, e.g. ACCESS1-0 (not the member's id)
            if !occursin(name, filename) && !(name in keys(MODEL_NAME_FIXES))
                @warn "model name as read from metadata of stored .nc file ($name) and used as dimension name is not identical to name appearing in its path ($filename)"
            end
            member_id_temp = memberIDsFromPaths([path])[1]
            model_id_temp = String(split(member_id_temp, MODEL_MEMBER_DELIM)[1])
            props["model_id"] = fixModelNameMetadata(model_id_temp)
            props["member_id"] = uniqueMemberID(attributes, props["model_id"])
            #props["physics"] = physicsFromMember(props["member_id"])
        elseif !is_model_data && dtype in [OBS_DATA, MODEL_OBS_DATA]
            props["model_id"] = filename
        else
            continue
        end
        # meta data that is stored in YAXArrays
        for key in collect(keys(dsVar.attrib)) #["units", "mip_era", "grid_label"]
            if haskey(attributes, key)
                props[key] = attributes[key]
            end
        end
        dimension_names = dimnames(dsVar)
        dimensions = Vector(undef, length(dimension_names))
        for (idx_dim, d) in enumerate(dimension_names)
            if d in ["bnds", "string21"]
                continue
            end
            if d == "time"
                # NOTE: just YEAR and MONTH are saved in the time dimension
                times = map(x -> DateTime(Dates.year(x), Dates.month(x)), dsVar[d][:])
                dimensions[idx_dim] = Dim{Symbol(d)}(collect(times))
            else
                dimensions[idx_dim] = Dim{Symbol(d)}(collect(dsVar[d][:]))
            end
        end
        push!(data, YAXArray(Tuple(dimensions), Array(dsVar), props))
        # replace missing values by NaN?
        #push!(data, YAXArray(Tuple(dimensions), coalesce.(Array(dsVar), NaN), props)
        close(ds)
    end
    return isempty(data) ? nothing : mergeDataFromMultipleFiles(data, sorted; meta=meta_info)
end



"""
    loadClimateData(
        all_paths::Vector{Vector{String}},
        ids::Vector{String};
        meta_data::Union{Vector{Dict{String, T}}, Nothing} = nothing,
        sorted::Bool = true,
        dtype::DataType = MODEL_OBS_DATA,
        fn_format::Symbol = :cmip
    ) where T <: Any

Load data at `all_paths`, every subvector refers the the paths of the data of a single 
dataset. 

Return a ClimateData instance when `dtype=MODEL_OBS_DATA` (default) and a DataMap with 
model (obs) data when `dtype=MODEL_DATA` (`dtype=OBS_DATA`).
"""
function loadClimateData(
    all_paths::Vector{Vector{String}},
    ids::Vector{String};
    meta_data::Union{Vector{Dict{String, T}}, Nothing} = nothing,
    sorted::Bool = true,
    dtype::DataType = MODEL_OBS_DATA,
    fn_format::Symbol = :cmip
) where T <: Any
    absent(x::Union{Vector, Nothing}) = isnothing(x) || isempty(x)
    if !absent(meta_data) && (length(all_paths) != length(meta_data))
        throw(ArgumentError("size of paths vector and meta data must be equal. Found: paths: $(length(all_paths)), meta_data: $(length(meta_data))"))
    end
    if length(all_paths) != length(ids)
        throw(ArgumentError("size of paths vector and ids must be equal. Found: paths: $(length(all_paths)), ids: $(length(ids))"))
    end
    getData(dtype::DataType) = map(
        (paths, meta) -> loadPreprocData(paths, fn_format; sorted, dtype, meta_info=meta),
        all_paths, 
        absent(meta_data) ? fill(nothing, length(all_paths)) : meta_data
    )
    models = dtype != OBS_DATA ? filter(!isnothing, getData(MODEL_DATA)) : []
    observations = dtype != MODEL_DATA ? filter(!isnothing, getData(OBS_DATA)) : []
    model_map = !isempty(models) ? buildDataMap(models, ids) : DataMap()
    obs_map = !isempty(observations) ? buildDataMap(observations, ids) : DataMap()
    return dtype == MODEL_DATA ? 
        model_map : 
        (dtype == OBS_DATA ? obs_map : ClimateData(model_map, obs_map))
end



"""
    loadDataFromESMValToolRecipes(
        path_data::String,
        path_recipes::String;
        dir_per_var::Bool = true,
        subset::Union{Dict, Nothing} = nothing,
        preview::Bool = false,
        sorted::Bool = true
    )

Load the data from the ESMValTool recipes at `path_recipes` or, if `preview` is true, load
meta data only.

# Arguments:
- `path_data`: top level directory were data is stored.
- `path_recipes`: directory were ESMValTool recipes are stored; these must be the versions 
in the run folder generated by ESMValTool named  RECIPENAME_filled.yml.
- `dir_per_var`: set to true (default) if there is a subdirectory in `path_data` for every 
climate variable to be loaded.
- `subset`: TODO!
- `preview`: set to true if only meta data should be loaded (default: false).
- `sorted`: if true (default), the data is sorted alphabetically wrt model names.
"""
function loadDataFromESMValToolRecipes(
    path_data::String,
    path_recipes::String;
    dir_per_var::Bool = true,
    subset::Union{Dict, Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true, 
    dtype::DataType = MODEL_OBS_DATA,
    fn_format::Symbol = :esmvaltool
)
    checkDataStructure(path_data, dir_per_var)
    meta_data = metaDataFromESMValToolRecipes(path_recipes; constraint=subset)
    paths = resolvePathsFromMetaData.(meta_data, path_data, dir_per_var; constraint=subset)
    if !isnothing(subset) && !isnothing(get(subset, "subset_shared", nothing))        
        paths = filterPathsSharedModels(paths, subset["subset_shared"])
    end
    if preview
        return collect(zip(meta_data, paths))
    end
    ids = map(x -> x.id, meta_data)
    return loadClimateData(
        paths, ids; sorted, dtype, fn_format, meta_data=metadataToDict.(meta_data)
    )
end


"""
    loadDataFromYAML(
        content::Dict;
        arg_constraint::Union{Dict, Nothing} = nothing,
        preview::Bool = false,
        sorted::Bool = true,
        dtype::DataType = MODEL_OBS_DATA
    )

Load a `ClimateData`-instance that contains the data specified in `content`, potentially 
constraint by values in `arg_constraint`.

# Arguments:
- `preview::Bool`: if true (default: false), return metadata and corresponding paths without 
actually loading any data.
"""
function loadDataFromYAML(
    content::Dict;
    arg_constraint::Union{Dict, Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true, 
    dtype::DataType = MODEL_OBS_DATA,
    fn_format::Symbol = :esmvaltool
)
    fn_err(x) = throw(ArgumentError("$(x) must be provided in config yaml file!"))
    datasets = get(() -> fn_err("datasets"), content, "datasets")
    base_path = get(() -> fn_err("base_path_data"), content, "base_path_data")

    all_meta = Vector{}(undef, length(datasets))
    all_paths = Vector{}(undef, length(datasets))
    for (i, ds) in enumerate(datasets)
        dir_per_var = get(ds, "dir_per_var", true)
        path_data = joinpath(base_path, get(() -> fn_err("base_dir"), ds, "base_dir"))
        checkDataStructure(path_data, dir_per_var)
        # merge arg_constraint has precedence over ds_subset here!
        ds_constraint = get(ds, "subset", nothing)
        if !isnothing(ds_constraint) && !isnothing(arg_constraint)
            warn_msg = "identical subset keys: arg constraint has precedence over constraint in yaml file!"
            ds_constraint = joinDicts(ds_constraint, arg_constraint; warn_msg) 
        end
        meta_data = metaDataFromYAML(ds)
        paths = resolvePathsFromMetaData.(meta_data, path_data, dir_per_var; constraint=ds_constraint)
        if !isnothing(ds_constraint) && !isnothing(get(ds_constraint, "subset_shared", nothing))        
            paths = filterPathsSharedModels(paths, ds_constraint["subset_shared"])
        end
        all_paths[i] = paths 
        all_meta[i] = meta_data 
    end
    paths = vcat(all_paths...)
    meta_data = vcat(all_meta...)
    # if argument provided, apply subset shared across all datasets
    level_all_data = get(arg_constraint, "subset_shared", nothing)
    if !isnothing(level_all_data)
        paths = filterPathsSharedModels(paths, level_all_data)
    end
    if preview
        return collect(zip(meta_data, paths))
    end
    ids = map(x -> x.id, meta_data)
    return loadClimateData(
        paths, ids; sorted, dtype, fn_format, meta_data=metadataToDict.(meta_data)
    )    
end

function loadDataFromYAML(
    path_config::String;
    arg_constraint::Union{Dict,Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true
)
    return loadDataFromYAML(YAML.load_file(path_config); arg_constraint, preview, sorted)
end

# Loading data directly from given directories
"""
    loadData(
        path_data_dir::String, 
        id::String;
        meta_data::Dict{String, T} = Dict{String, Any}(),
        constraint::Union{Dict, Nothing} = nothing,
        sorted::Bool = true, 
        dtype::DataType = MODEL_OBS_DATA,
        fn_format::Symbol=:cmip
    ) where T <: Any

Return DataMap with entry `id` with the data (model and/or obs depending on `dtype`) from 
all .nc files in `paths` and all .nc files in all directories in `paths`, possibly 
constraint by `constraint`.
"""
function loadData(
    paths::Vector{String}, 
    id::String;
    meta_data::Dict{String, T} = Dict{String, Any}(),
    constraint::Union{Dict, Nothing} = nothing,
    sorted::Bool = true, 
    dtype::DataType = MODEL_OBS_DATA,
    fn_format::Symbol=:cmip
) where T <: Any
    paths_dirs = filter(x -> isdir(x), paths)
    paths_ncfiles = filter(x -> isfile(x) && endswith(x, ".nc"), paths)
    paths_to_files = vcat(collectNCFilePaths.(paths_dirs; constraint)..., paths_ncfiles)
    return loadClimateData([paths_to_files], [id]; meta_data = [meta_data], sorted, dtype, fn_format)    
end


"""
    loadData(
        paths_data_dirs::Vector{String}, 
        data_ids::Vector{String};
        meta_data::Vector{Dict{String, T}} = Vector{Dict{String, Any}}(),
        constraint::Union{Dict, Nothing} = nothing,
        sorted::Bool = true, 
        dtype::DataType = MODEL_OBS_DATA,
        fn_format::Symbol=:cmip
    ) where T <: Any

Return DataMap with entries `data_ids` with the data (model and/or obs depending on `dtype`) 
from all .nc files in all directories in `paths_dirs`, possibly constraint by `constraint`.
"""
function loadData(
    paths_data_dirs::Vector{String}, 
    data_ids::Vector{String};
    meta_data::Vector{Dict{String, T}} = Vector{Dict{String, Any}}(),
    constraint::Union{Dict, Nothing} = nothing,
    sorted::Bool = true, 
    dtype::DataType = MODEL_OBS_DATA,
    fn_format::Symbol=:cmip
) where T <: Any
    paths_to_files = collectNCFilePaths.(paths_data_dirs; constraint)
    if !isnothing(constraint) && !isempty(constraint) && !isnothing(get(constraint, "subset_shared", nothing))        
        paths_to_files = filterPathsSharedModels(paths_to_files, constraint["subset_shared"])
    end
    return loadClimateData(paths_to_files, data_ids; sorted, dtype, fn_format, meta_data)
end


"""
    loadData(
        paths_data_dirs::Vector{String}, 
        data_ids::Vector{String};
        meta_data::Vector{Dict{String, T}} = Vector{Dict{String, Any}}(),
        constraint::Union{Dict, Nothing} = nothing,
        sorted::Bool = true, 
        dtype::DataType = MODEL_OBS_DATA,
        fn_format::Symbol=:cmip
    ) where T <: Any

Return DataMap with entries `data_ids` with the data (model and/or obs depending on `dtype`) 
from all .nc files in all directories in `paths_dirs`, possibly constraint by `constraint`.

For every loaded dataset (entry in built DataMap), files are loaded from several directories;
each entry in `paths_data_dirs` points the the vector of data directories from where data 
is loaded for that dataset.
"""
function loadData(
    paths_data_dirs::Vector{Vector{String}}, 
    data_ids::Vector{String};
    meta_data::Vector{Dict{String, T}} = Vector{Dict{String, Any}}(),
    constraint::Union{Dict, Nothing} = nothing,
    sorted::Bool = true, 
    dtype::DataType = MODEL_OBS_DATA,
    fn_format::Symbol=:cmip
) where T <: Any
    paths_to_files = map(paths_data_dirs) do paths 
        vcat(collectNCFilePaths.(paths; constraint)...)
    end
    if !isnothing(constraint) && !isempty(constraint) && !isnothing(get(constraint, "subset_shared", nothing))        
        paths_to_files = filterPathsSharedModels(paths_to_files, constraint["subset_shared"])
    end
    return loadClimateData(paths_to_files, data_ids; sorted, dtype, fn_format, meta_data)
end
