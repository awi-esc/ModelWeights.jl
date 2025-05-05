using Dates
using DimensionalData
using NCDatasets
using Setfield
using YAML
using YAXArrays

"""
    getUniqueMemberIds(meta::Dict, model_names::Vector{String})

Create a vector that contains for every model in `model_names` the unique ids
of its members, each with the model name followed by '#' as prefix .

The unique member ids correspond to the variant labels of CMIP6 models, e.g. r1i1p1f1.

# Arguments:
- `meta::Dict`: For CMIP5 data, must have keys: 'mip_era', 'realization',
'initialization_method', 'physics_version'. For CMIP6 data must have keys: 
'variant_label', 'grid_label'.
- `model_names::Vector{String}`: Vector of strings containing model names for 
every model member, i.e. length is sum of the number of members over all models.
"""
function getUniqueMemberIds(meta::Dict, model_names::Vector{String})
    meta_subdict = Dict{String,Vector}()
    keys_model_ids = [
        "realization",
        "physics_version",
        "initialization_method",
        "mip_era",
        "grid_label",
        "variant_label",
    ]
    n_members = length(model_names)
    members = Vector{String}(undef, n_members)
    for key in filter(k -> k in keys_model_ids, keys(meta))
        val = meta[key]
        if isa(val, String)
            meta_subdict[key] = repeat([val], outer = n_members)
        else
            meta_subdict[key] = deepcopy(val)
        end
    end
    mip_eras = meta_subdict["mip_era"]
    indices_cmip5 = findall(x -> !ismissing(x) && x == "CMIP5", mip_eras)
    indices_cmip6 = findall(x -> !ismissing(x) && x == "CMIP6", mip_eras)
    if !isempty(indices_cmip5)
        variants = buildCMIP5EnsembleMember(
            meta_subdict["realization"][indices_cmip5],
            meta_subdict["initialization_method"][indices_cmip5],
            meta_subdict["physics_version"][indices_cmip5],
        )
        members[indices_cmip5] = map(
            x -> join(x, MODEL_MEMBER_DELIM, MODEL_MEMBER_DELIM),
            zip(model_names[indices_cmip5], variants),
        )
        @debug "For CMIP5, full model names dont include grid."
    end
    if !isempty(indices_cmip6)
        variants = meta_subdict["variant_label"][indices_cmip6]
        grids = meta_subdict["grid_label"][indices_cmip6]
        members[indices_cmip6] = map(
            x -> join(x, MODEL_MEMBER_DELIM, "_"),
            zip(model_names[indices_cmip6], variants, grids),
        )
    end
    # from vector of strings of length of sum over members of all models make
    # a vector of vectors where each subvector contains the unique ids of the
    # respective model's members
    models_uniq = unique(model_names)
    n_models = length(models_uniq)
    result = Vector{Vector{String}}(undef, n_models)
    for (i, model) in enumerate(models_uniq)
        indices = findall(x -> x == model, model_names)
        result[i] = members[indices]
    end
    return vcat(result...)
end


"""
    filterPaths(paths::Vector{String}, shared_models::Vector{String})

Return new vector with paths in `paths` that point to models in `shared_models`.

# Arguments:
- `shared_models`: have form 'modelname#memberID[_grid]' on level of model
members or just 'modelname' on level of models.
"""
function filterPaths(paths::Vector{String}, shared_models::Vector{String})
    subset_paths = Vector{String}()
    for path in paths
        keep = applyModelConstraints(path, shared_models)
        if keep
            push!(subset_paths, path)
        end
    end
    return subset_paths
end


"""
    subsetModelData(data::AbstractArray, shared_models::Vector{String})

Return a YAXArray that is the subset of `data` containing only data from 
the models specified in `shared_models`. 

Takes care of metadata.

# Arguments:
- `data`: must have dimension 'member' or 'model'
- `shared_models`: vector of models, which can either be on level of models 
or members of models ('modelname#memberID[_grid]').
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
            models_data = getModelsFromMemberIDs(dim_names) # NOTE: should yield same: models_data = data.properties["model_names"]
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
        #println("$(key): $(length(data.properties[key]))")
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
    shared_models = level == MEMBER ? getSharedMembers(datamap) : getSharedModels(datamap)
    if isempty(shared_models)
        @warn "Vector of models to subset data to is empty!"
        return datamap
    end
    subset = DataMap()
    for data in all_data
        df = deepcopy(data)
        subset[df.properties["_id"]] = subsetModelData(df, shared_models)
    end
    return subset
end


"""
    loadPreprocData(meta::Dict{String, Any}, is_model_data::Bool)

Create a tuple with a vector of YAXArrays for data specified in `meta` and a 
Dictionary containing the metadata of all loaded data. 

Load the data from paths specified in `meta.paths` and create a meta dictionary
that contains the metadata keys from every loaded dataset. Each key maps to a vector 
of values, one for each loaded dataset, which is set to missing if that key 
hadn't been present in this datasets own metadata.
"""
function loadPreprocData(meta::Dict{String,Any}, is_model_data::Bool)
    n_files = length(meta["_paths"])
    data = Vector{YAXArray}(undef, n_files)
    meta_dict = Dict{String,Any}()
    source_names = repeat([""], outer = n_files)

    filenames = first.(splitext.(basename.(meta["_paths"])))
    climVars = first.(split.(basename.(dirname.(meta["_paths"])), "_"))
    for (i, file) in enumerate(meta["_paths"])
        @debug "processing file.." * file
        filename = filenames[i]
        climVar = climVars[i]

        ds = NCDataset(file)
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
        if is_model_data
            # add mip_era for models since it is not provided in CMIP5-models
            model_key = getCMIPModelsKey(Dict(ds.attrib))
            get!(attributes, "mip_era", "CMIP5")
            # check model names as retrieved from the metadata for potential inconsistencies wrt filename
            name = ds.attrib[model_key] # just the model name, e.g. ACCESS1-0 (not the member's id)
            if !occursin(name, filename) && !(name in keys(MODEL_NAME_FIXES))
                @warn "model name as read from metadata of stored .nc file ($name) and used as dimension name is not identical to name appearing in its path ($filename)"
            end
            model_id = getMemberIDsFromPaths([file])[1]
            source_names[i] = split(model_id, MODEL_MEMBER_DELIM)[1]
        else
            source_names[i] = filename
        end
        # update metadata-dictionary for all processed files with the
        # metadata from the current file

        keys_remain = filter(x -> x in KEYS_METADATA, keys(attributes))
        for key in keys_remain #keys(attributes)
            values =
                get!(meta_dict, key, repeat(Union{Missing,Any}[missing], outer = n_files))
            values[i] = attributes[key]
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
                # times = 1:length(dsVar[d][:])
                dimensions[idx_dim] = Dim{Symbol(d)}(collect(times))
            else
                dimensions[idx_dim] = Dim{Symbol(d)}(collect(dsVar[d][:]))
            end
        end
        data[i] = YAXArray(Tuple(dimensions), Array(dsVar))
        # replace missing values by NaN?
        #data[i] = YAXArray(Tuple(dimensions), coalesce.(Array(dsVar), NaN))
        close(ds)
    end
    updateMetadata!(meta_dict, collect(source_names), is_model_data)
    return (data, meta_dict)
end


"""
    mergeLoadedData(
        data_vec::Vector{YAXArray}, 
        meta_dict::Dict{String, Any}, 
        is_model_data::Bool
    )

All data is assumed to be defined on the same grid.
"""
function mergeLoadedData(
    data_vec::Vector{YAXArray},
    meta_dict::Dict{String,Any},
    is_model_data::Bool,
)
    if length(data_vec) == 0
        @warn "No data loaded!"
        return nothing
    end
    #dimData = cat(data_vec..., dims=3) # way too slow!
    data_sizes = unique(map(size, data_vec))
    if length(data_sizes) != 1
        if !all(map(x -> hasdim(x, :time), data_vec))
            msg = "Data does not have the same size across all models: $(data_sizes)"
            throw(ArgumentError(msg))
        else
            # if difference only in time, use maximal possible timeseries and add NaNs
            year_min = minimum(map(x -> minimum(map(Dates.year, dims(x, :time))), data_vec))
            year_max = maximum(map(x -> maximum(map(Dates.year, dims(x, :time))), data_vec))
            nb_years = year_max - year_min + 1

            timerange = DateTime(year_min):Year(1):DateTime(year_max)
            for (i, ds) in enumerate(data_vec)
                s = map(length, otherdims(ds, :time))
                dat = Array{eltype(ds)}(undef, s..., nb_years) # if ds allows missing values, undef is initialized with missing
                ds_extended =
                    YAXArray((otherdims(ds, :time)..., Dim{:time}(timerange)), dat)
                ds_extended[time=Where(
                    x->Dates.year(x) in map(Dates.year, dims(ds, :time)),
                )] = ds # ds[time=:]
                data_vec[i] = ds_extended
            end
        end
    end
    var_axis =
        is_model_data ? Dim{:member}(meta_dict["member_names"]) :
        Dim{:model}(meta_dict["model_names"])
    dimData = concatenatecubes(data_vec, var_axis)
    dimData = YAXArray(dimData.axes, dimData.data, meta_dict)

    # Sanity checks that no dataset exists more than once
    if is_model_data
        members = dims(dimData, :member)
        if length(members) != length(unique(members))
            duplicates = [m for m in members if sum(members .== m) > 1]
            #paths = [filter(x -> occursin(split(dup, MODEL_MEMBER_DELIM)[1], x), meta.paths) for dup in duplicates]
            @warn "Some datasets appear more than once" duplicates
        end
    end
    return dimData
end


"""
    loadDataFromMetadata(
        meta_data::Dict{String, Dict{String, Any}}, is_model_data::Bool
    )

"""
function loadDataFromMetadata(meta_data::Dict{String,Dict{String,Any}}, is_model_data::Bool)
    results = DataMap()
    for (id, meta) in meta_data
        # loads data at level of model members
        @info "load $id"
        data_vec, meta_dict = loadPreprocData(meta, is_model_data)
        data = mergeLoadedData(data_vec, meta_dict, is_model_data)
        merge!(data.properties, meta)
        results[meta["_id"]] = data
    end
    return results
end


"""
    getMemberIDsFromPaths(all_paths::Vector{String})

For every path in `all_paths` return a string of the form modelname#memberID[_grid]
that identifies the corresponding model member.

The abbreviation of the grid is added to the model name for CMIP6 models.
"""
function getMemberIDsFromPaths(all_paths::Vector{String})
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
    searchModelInPaths(model::String, paths::Vector{String})

# Arguments:
- `model_id::String`: has form modelname#memberID[_grid]
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
    for (_, meta_path) in enumerate(paths)
        found_model = occursin("_" * model * "_", meta_path)
        if !found_model
            continue
        else
            if has_member
                found_member = occursin("_" * member * "_", meta_path)
                if found_member
                    if has_grid
                        found_grid = occursin("_" * grid, meta_path)
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
    getSharedModelsFromPaths(
        all_paths::Vector{Vector{String}}, all_models::Vector{String}
    )
     
Return vector with those models in 'all_models' for which a path is given in 
every subvector of 'meta_data'.
"""
function getSharedModelsFromPaths(
    all_paths::Vector{Vector{String}},
    all_models::Vector{String},
)
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
            @debug "Dictionary contains keys 'source_id' used in CMIP6 ($(meta["source_id"])) and 'model_id' used in CMIP5 ($(meta["model_id"])). 'source_id' is used!"
        end
        return "source_id"
    elseif "model_id" in attributes
        return "model_id"
    else
        msg = "Metadata must contain one of 'source_id' (pointing to names of CMIP6 models or 'model_id' used for CMIP5 models."
        throw(ArgumentError(msg))
    end
end



"""
    filterPathsSharedModels!(
        meta_data::Dict{String, Dict{String, Any}},
        subset_shared::Union{LEVEL, Nothing}
    )

    In particular important before loading data if data should be subset
"""
function filterPathsSharedModels!(
    meta_data::Dict{String,Dict{String,Any}},
    subset_shared::Union{LEVEL,Nothing},
)
    all_paths = map(x -> x["_paths"], values(meta_data))
    all_models = getMemberIDsFromPaths(vcat(all_paths...))
    if subset_shared == MODEL
        all_models = unique(getModelsFromMemberIDs(all_models))
    end
    shared_models = getSharedModelsFromPaths(all_paths, all_models)
    if isempty(shared_models)
        @warn "No models shared across data!"
    end
    for (id, meta) in meta_data
        meta_data[id]["_paths"] = filterPaths(meta["_paths"], shared_models)
    end
    return nothing
end

function getModelsFromMemberIDs(members::Vector{String})
    return map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), members)
end


function getPhysicsFromMembers(members::Vector{String})
    regex = r"(p\d+)(f\d+)?(_.*)?$"
    return String.(map(x -> match(regex, x).captures[1], members))
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
    models = unique(getModelsFromMemberIDs(members))
    for model in models
        # retrieve allowed physics as given in members 
        member_ids = filter(m -> startswith(m, model * MODEL_MEMBER_DELIM), members)
        physics = getPhysicsFromMembers(member_ids)
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
                data[ds.properties["_id"]] = subsetModelData(ds, members_kept)
            end
        end
    end
    if !isnothing(subset_shared)
        shared_models =
            subset_shared == MEMBER ? getSharedMembers(data) : getSharedModels(data)
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

For each model and variable (if several given), compute a summary statistic 
(default: mean) across all members of that model. Instead of 'member', the 
returned YAXArray has dimension 'model'.

# Arguments:
- `data::YAXArray`: a YAXArray with at least dimension 'member'
- `updateMeta::Bool`: set true if the vectors in the metadata refer to 
different models. Set to false if vectors refer to different variables.
"""
function summarizeEnsembleMembersVector(
    data::YAXArray,
    updateMeta::Bool;
    fn::Function = Statistics.mean,
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
        dat = data[model=Where(x->x==m)]
        average =
            isempty(dimensions) ? fn(dat) : mapslices(x -> fn(x), dat; dims = (:model,))
        summarized_data[model=At(m)] = average
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
