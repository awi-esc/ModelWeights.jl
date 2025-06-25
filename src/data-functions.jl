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
    return sharedModelsFromPaths(all_paths, all_models)
end


"""
    loadPreprocData(
        paths::Vector{String}; 
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
    paths::Vector{String}; 
    sorted::Bool = true, 
    dtype::DataType = MODEL_OBS_DATA, 
    meta::Union{Dict{String, T}, Nothing} = nothing
) where T <: Any
    data = Vector{YAXArray}()    
    filenames = first.(splitext.(basename.(paths)))
    climVars = first.(split.(basename.(dirname.(paths)), "_"))
    for (i, path) in enumerate(paths)
        @debug "processing file.." * path
        filename = filenames[i]
        climVar = climVars[i]
        is_model_data = any(map(x -> startswith(lowercase(filename), x), ["cmip5", "cmip6"]))

        ds = NCDataset(path)
        dsVar = ds[climVar]
        # if climVar == "amoc"
        #     dsVar = ds["msftmz"]
        # end
        attributes = merge(dsVar.attrib, ds.attrib)
        if climVar == "msftmz"
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
    return isempty(data) ? nothing : mergeDataFromMultipleFiles(data, sorted; meta)
end


function loadClimateData(
    all_paths::Vector{Vector{String}},
    ids::Vector{String};
    meta_data::Union{Vector{Dict{String, T}}, Nothing} = nothing,
    sorted::Bool = true,
    dtype::DataType = MODEL_OBS_DATA
) where T <: Any
    absent(x::Union{Vector, Nothing}) = isnothing(x) || isempty(x)
    if !absent(meta_data) && (length(all_paths) != length(meta_data))
        throw(ArgumentError("size of paths vector and meta data must be equal. Found: paths: $(length(all_paths)), meta_data: $(length(meta_data))"))
    end

    getData(dtype::DataType) = map(
        (paths, meta) -> loadPreprocData(paths; sorted, meta, dtype),
        all_paths, 
        absent(meta_data) ? fill(nothing, length(all_paths)) : meta_data
    )
    models = dtype != OBS_DATA ? filter(!isnothing, getData(MODEL_DATA)) : []
    observations = dtype != MODEL_DATA ? filter(!isnothing, getData(OBS_DATA)) : []
    model_map = !isempty(models) ? buildDataMap(models, ids) : DataMap()
    obs_map = !isempty(observations) ? buildDataMap(observations, ids) : DataMap()
    return ClimateData(model_map, obs_map)
end


function filterPathsSharedModels!(all_paths::Vector{Vector{String}}, subset_shared::Level)
    filterPathsSharedModels!(all_paths, sharedModels(all_paths, subset_shared))
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
    dtype::DataType = MODEL_OBS_DATA
)
    checkDataStructure(path_data, dir_per_var)
    meta_data = metaDataFromESMValToolRecipes(path_recipes; constraint = subset)   
    paths = resolvePathsFromMetaData(meta_data, path_data; dir_per_var, subset)
    if preview
        return collect(zip(meta_data, paths))
    end
    ids = map(x -> x.id, meta_data)
    return loadClimateData(paths, ids; sorted, dtype, meta_data=metadataToDict.(meta_data))
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
    dtype::DataType = MODEL_OBS_DATA
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
        paths = resolvePathsFromMetaData(meta_data, path_data; dir_per_var, subset=ds_constraint)
        all_paths[i] = paths 
        all_meta[i] = meta_data 
    end
    # apply subset shared across all datasets
    level_all_data = get(arg_constraint, "subset_shared", nothing)
    if !isnothing(level_all_data)
        shared_models = sharedModels(vcat(all_paths...), level_all_data)
        map(all_paths) do paths 
            filterPathsSharedModels!(paths, shared_models)
        end
    end
    paths = vcat(all_paths...)
    meta_data = vcat(all_meta...)
    if preview
        return collect(zip(meta_data, paths))
    end
    ids = map(x -> x.id, meta_data)
    return loadClimateData(paths, ids; sorted, dtype, meta_data=metadataToDict.(meta_data))    
end

function loadDataFromYAML(
    path_config::String;
    arg_constraint::Union{Dict,Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true
)
    return loadDataFromYAML(YAML.load_file(path_config); arg_constraint, preview, sorted)
end

# Loading data for given paths
function loadDataFromDirs(
    paths_data_dirs::Vector{Vector{String}}, 
    data_ids::Vector{String};
    meta_dicts::Vector{Dict{String, T}} = Vector{Dict{String, Any}}(),
    constraint::Union{Dict, Nothing} = nothing,
    sorted::Bool = true, 
    dtype::DataType = MODEL_OBS_DATA
) where T <: Any
    paths = map(paths_data_dirs) do paths 
        vcat(collectNCFilePaths.(paths; constraint)...)
    end
    return loadClimateData(paths, data_ids; meta_data=meta_dicts, sorted, dtype)
end

