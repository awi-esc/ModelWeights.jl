"""
    uniqueMemberID(meta::Dict, model::String)

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
function uniqueMemberID(meta::Dict, model::String)
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
        return nothing
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
    subsetModelData(datamap::DataMap; level::LEVEL=MEMBER)

For those datasets in `datamap` that specify data on the level `level` 
(i.e. have dimension :member or :model), return a new DataMap with subset of 
data s.t. the new datasets all have the same models (level=MODEL) or members 
(level=MEMBER).

If no models are shared across datasets, return the input `datamap`.
"""
function subsetModelData(datamap::DataMap; level::LEVEL = MEMBER)
    all_data = collect(values(datamap))
    shared_models = level == MEMBER ? sharedMembers(datamap) : sharedModels(datamap)
    if isempty(shared_models)
        @warn "Vector of models to subset data to is empty!"
        return datamap
    end
    subset = DataMap()
    for data in all_data
        df = deepcopy(data)
        subset[df.properties["id"]] = subsetModelData(df, shared_models)
    end
    return subset
end


"""
    loadPreprocData(meta::MetaData; sorted::Bool = true)

Return data loaded from paths specified in `meta.paths` as single YAXArray.  The names of
the data files are assumed to follow the CMIP-standard: TODO!
"""
function loadPreprocData(meta::MetaData; sorted::Bool = true)
    n_files = length(meta.paths)
    data = Vector{YAXArray}(undef, n_files)
    
    filenames = first.(splitext.(basename.(meta.paths)))
    climVars = first.(split.(basename.(dirname.(meta.paths)), "_"))
    for (i, path) in enumerate(meta.paths)
        @debug "processing file.." * path
        filename = filenames[i]
        climVar = climVars[i]

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
        props = Dict{String, Union{String, Number, Missing}}() # just for this file!
        props["path"] = path

        if meta.is_model_data
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
        else
            props["model_id"] = filename
        end

        for key in ["units", "mip_era", "grid_label"]
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
        data[i] = YAXArray(Tuple(dimensions), Array(dsVar), props)
        # replace missing values by NaN?
        #data[i] = YAXArray(Tuple(dimensions), coalesce.(Array(dsVar), NaN))
        close(ds)
    end
    return mergeDataFromMultipleFiles(data, meta, sorted)
end


"""
    mergeDataFromMultipleFiles(
        data::Vector{YAXArray}, meta_dict::Dict{String, Vector}, meta::MetaData, sorted::Bool
    )

Combine arrays in `data` into a single YAXArray with meta combined from all datasets into lists, 
with missing if key wasnt in a dataset.
"""
function mergeMetaDataFromMultipleFiles(data::Vector{YAXArray})
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
    mergeDataFromMultipleFiles(data::Vector{YAXArray}, meta::MetaData, sorted::Bool)

Combine arrays in `data` into a single YAXArray with meta combined from all datasets into lists, 
with missing if key wasnt in a dataset, and with properties in `meta` (except for paths). 
If `sorted` is true, model dimension of returned data is sorted alphabetically and the 
entries in the metadata dictionary of the returned array which are vectors are sorted accordingly.

All data must be defined on the same grid. For timeseries data, the time dimension might cover 
different ranges, in that case they maximal timeseries is used and filled with NaN for missing values. 
"""
function mergeDataFromMultipleFiles(data::Vector{YAXArray}, meta::MetaData, sorted::Bool)
    if isempty(data)
        @warn "No data loaded!"
        return nothing
    end
    if !all(x -> haskey(x.properties, "model_id"), data)
        throw(ArgumentError("All data must have model_id in metadata!"))
    end
    #dimData = cat(data..., dims=3) # way too slow!
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

    hasMembers = all(x -> haskey(x.properties, "member_id"), data)
    dim = hasMembers ? :member : :model
    names = map(x -> x.properties[String(dim) * "_id"], data)
    if sorted
        sort_indices = sortperm(names)
        names = names[sort_indices]
        data = data[sort_indices]
    end
    meta_dict = mergeMetaDataFromMultipleFiles(data)    
    merge!(meta_dict, metadataToDict(meta; exclude=[:paths]))
    if hasMembers
        renameDictKeys!(
            meta_dict, 
            [("model_id", "model_names"), ("member_id", "member_names"), ("path", "paths")]
        )
    end
    dimData = concatenatecubes(data, Dim{dim}(names))
    dimData = YAXArray(dimData.axes, dimData.data, meta_dict)
    # Sanity checks that no dataset exists more than once
    if hasMembers
        members = dims(dimData, :member)
        if length(members) != length(unique(members))
            duplicates = unique([m for m in members if sum(members .== m) > 1])
            @warn "Some datasets appear more than once" duplicates
        end
    end
    return dimData
end


"""
    loadDataFromMetadata(meta_data::Dict{String, MetaData}, sorted::Bool)

Load a `DataMap`-instance that contains the data on level of model members specified in `meta_data`.

# Arguments:
- `meta_data::Dict`: keys are ids of datasets to be loaded, values is their MetaData
- `sorted::Bool`: set to true if data should be sorted alphabetically by model names.
"""
function loadDataFromMetadata(meta_data::MetaDataMap, sorted::Bool)
    @debug "retrieved metadata:" meta_data
    if isempty(meta_data)
        throw(ArgumentError("No metadata provided!"))
    end
    results = DataMap()
    for (id, meta) in meta_data
        @info "load Data for: $id"
        data = loadPreprocData(meta; sorted)
        results[meta.id] = data
    end
    return results
end


"""
    memberIDsFromPaths(all_paths::Vector{String})

For every path in `all_paths` return a string of the form modelname#memberID[_grid]
that identifies the corresponding model member.

The abbreviation of the grid is added to the model name for CMIP6 models.
"""
function memberIDsFromPaths(all_paths::Vector{String})
    all_filenames = split.(basename.(all_paths), "_")
    # model names are at predefined position in filenames (ERA_name_mip_exp_id_variable[_grid].nc)
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
    searchModelInPaths(model_id::String, paths::Vector{String})

Return true if `model_id` is found in filename of any path in `paths`, else false.

The paths are assumed to follow the standard CMIP filename structure, i.e. <variable_id>_<table_id>_<source_id>_<experiment_id >_<member_id>_<grid_label>[_<time_range>].nc([see here](https://docs.google.com/document/d/1h0r8RZr_f3-8egBMMh7aqLwy3snpD6_MrDz1q8n5XUk/edit?tab=t.0)).

# Arguments:
- `model_id::String`: has form modelname[#memberID[_grid]]
- `paths::Vector{String}`: paths to be searched
"""
function searchModelInPaths(model_id::String, paths::Vector{String})
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
    sharedModelsFromPaths(all_paths::Vector{Vector{String}}, all_models::Vector{String})
     
Return vector with models in `all_models` for which a path is given in every subvector of `all_paths`.
"""
function sharedModelsFromPaths(all_paths::Vector{Vector{String}}, all_models::Vector{String})
    indices_shared = []
    for (idx, model) in enumerate(all_models)
        is_found = false
        for paths in all_paths
            is_found = searchModelInPaths(model, paths)
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
    filterPathsSharedModels!(
        meta_data::Dict{String, Dict{String, T}}, subset_shared::LEVEL
    ) where T <: Any

For every dataset in `meta_data`, filter '_paths' s.t. the paths vector for each dataset 
only contains models or model members (given by `subset_shared`) that are shared across all 
datasets.

# Arguments:
- `meta_data`: contains a metadata dictionary for every dataset.
- `subset_shared`: level on which datasets must match: they include data from the  same 
models (MODEL) or from the same model members (MEMBER).
"""
function filterPathsSharedModels!(meta_data::MetaDataMap, subset_shared::LEVEL)
    all_paths = map(x -> x.paths, values(meta_data))
    all_models = memberIDsFromPaths(vcat(all_paths...))
    if subset_shared == MODEL
        all_models = unique(modelsFromMemberIDs(all_models))
    end
    shared_models = sharedModelsFromPaths(all_paths, all_models)
    if isempty(shared_models)
        @warn "No models shared across data!"
    end
    for (id, meta) in meta_data
        mask = map(p -> applyModelConstraints(p, shared_models), meta.paths)
        meta_data[id].paths = meta.paths[mask]
    end
    return nothing
end

function filterPathsSharedModels!(meta_data::MetaDataMap, subset_shared::String)
    level = get(() -> throw(ArgumentError("subset_shared can only be one of: $(keys(LEVEL_LOOKUP))")), 
        LEVEL_LOOKUP, lowercase(subset_shared))
    filterPathsSharedModels!(meta_data, level)
    return nothing
end


function modelsFromMemberIDs(members::AbstractVector{<:String}; uniq::Bool=false)
    models = map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), members)
    return uniq ? unique(models) : models
end


function physicsFromMember(member::String)
    regex = r"(p\d+)(f\d+)?(_.*)?$"
    return String.(match(regex, member).captures[1])
end


"""
    alignPhysics(
        data::DataMap,
        members::Vector{String}, 
        subset_shared::Union{LEVEL, Nothing} = nothing)
    )

Return new DataMap with only the models retained that share the same physics as 
the respective model's members in `members`.

If `subset_shared` is set, resulting DataMap is subset accordingly.
"""
function alignPhysics(
    datamap::DataMap,
    members::Vector{String};
    subset_shared::Union{LEVEL,Nothing} = nothing,
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
        shared_models =
            subset_shared == MEMBER ? sharedMembers(data) : sharedModels(data)
        for (id, model_data) in data
            data[id] = subsetModelData(model_data, shared_models)
        end
    end
    return data
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
    setToSummarizedMembers!(data::DataMap)

Set values for every dataset in `data` to the average across all members of 
each model.
"""
function setToSummarizedMembers!(data::DataMap)
    for (k, current_data) in data
        if hasdim(current_data, :member)
            @info "average ensemble members for $k"
            data[k] = summarizeEnsembleMembersVector(current_data, true)
        end
    end
    return nothing
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


function addMasks!(datamap::DataMap, id_orog_data::String)
    orog_data = datamap[id_orog_data]
    datamap["mask_land"] = getMask(orog_data; mask_out_land = false)
    datamap["mask_ocean"] = getMask(orog_data; mask_out_land = true)
    return nothing
end


function checkDataStructure(path_data::String, dir_per_var::Bool, is_model_data::Bool)
    target_folder = basename(path_data)
    if dir_per_var && target_folder == "preproc"
        throw(ArgumentError("If dir_per_var is true, path_data should point to a parentdirectory of directory with name 'preproc'"))
    elseif !dir_per_var && is_model_data && target_folder != "preproc"
        throw(ArgumentError("If dir_per_var is false, path_data should point to folder with name 'preproc'"))
    end
    return nothing
end


function loadMetaData(
    meta_attribs::Vector{MetaData},
    path_data::String;
    dir_per_var::Bool = true,
    subset::Union{Dict, Nothing} = nothing
)
    if !isnothing(subset)
        constrainMetaData!(meta_attribs, subset)
    end
    addPaths!(meta_attribs, path_data, dir_per_var; constraint = subset)
    meta_data = metaVecToDict(meta_attribs)
    if !isnothing(subset) && !isnothing(get(subset, "subset_shared", nothing))        
        filterPathsSharedModels!(meta_data, subset["subset_shared"])
    end
    return meta_data
end


"""
    loadDataFromESMValToolRecipes(
        path_data::String,
        path_recipes::String;
        dir_per_var::Bool = true,
        is_model_data::Bool = true,
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
- `is_model_data`: set to true (default) for loading CMIP data, to false for observational 
data.
- `subset`: TODO!
- `preview`: set to true if only meta data should be loaded (default: false).
- `sorted`: if true (default), the data is sorted alphabetically wrt model names.
"""
function loadDataFromESMValToolRecipes(
    path_data::String,
    path_recipes::String;
    dir_per_var::Bool = true,
    is_model_data::Bool = true,
    subset::Union{Dict, Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true
)
    checkDataStructure(path_data, dir_per_var, is_model_data)
    attribs = metaAttributesFromESMValToolRecipes(path_recipes, is_model_data)
    meta_data = loadMetaData(attribs, path_data; dir_per_var, subset)
    return preview ? meta_data : loadDataFromMetadata(meta_data, sorted)
end


"""
    loadDataFromYAML(
        content::Dict;
        arg_constraint::Union{Dict, Nothing} = nothing,
        preview::Bool = false,
        sorted::Bool = true
    )

Load a `DataMap`-instance that contains the data specified in `content`, potentially 
constraint by values in `arg_constraint`.

# Arguments:
- `preview::Bool`: if true (default: false), return metadata without actually 
loading any data.
"""
function loadDataFromYAML(
    content::Dict;
    arg_constraint::Union{Dict, Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true
)
    fn_err(x) = throw(ArgumentError("$(x) must be provided in config yaml file!"))
    datasets = get(() -> fn_err("datasets"), content, "datasets")
    base_path = get(() -> fn_err("base_path_data"), content, "base_path_data")

    data = preview ? MetaDataMap() : DataMap()
    meta_data = MetaDataMap()
    for ds in datasets
        dir_per_var = get(ds, "dir_per_var", true)
        is_model_data = get(ds, "is_model_data", true)
        path_data = joinpath(base_path, get(() -> fn_err("base_dir"), ds, "base_dir"))
        checkDataStructure(path_data, dir_per_var, is_model_data)
        attribs = metaAttributesFromYAML(ds)
        # merge arg_constraint has precedence over ds_subset here!
        ds_constraint = get(ds, "subset", nothing)
        if !isnothing(ds_constraint) && !isnothing(arg_constraint)
            ds_constraint = joinDicts(ds_constraint, arg_constraint)
        end
        meta = loadMetaData(attribs, path_data; dir_per_var, subset = ds_constraint)
        meta_data = joinDataMaps(meta_data, meta)
    end
    
    level_all_data = get(arg_constraint, "subset_shared", nothing)
    if !isnothing(level_all_data)
        filterPathsSharedModels!(meta_data, level_all_data)
    end
    
    return preview ? meta_data : loadDataFromMetadata(meta_data, sorted)
end

function loadDataFromYAML(
    path_config::String;
    arg_constraint::Union{Dict,Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true
)
    return loadDataFromYAML(YAML.load_file(path_config); arg_constraint, preview, sorted)
end
