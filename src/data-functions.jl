using DimensionalData
using NCDatasets
using YAML

"""
    getUniqueMemberIds(meta::Dict, model_names::Vector{String})

Combine name of the model with the id of the respective model members. The
returned vector has length of the number of models and contains for every model
a vector with the unique ids of its members.

# Arguments:
- `meta`: for CMIP5 data, must have keys: 'mip_era', 'realization', 'initialization_method',
'physics_version'. For CMIP6 data must have keys: 'variant_label', 'grid_label'
- `model_names`: Vector of strings containing model names for every model member, i.e.
length is sum of the number of members over all models
"""
function getUniqueMemberIds(meta::Dict, model_names::Vector{String})
    meta_subdict = Dict{String, Vector}()
    keys_model_ids = [
        "realization", "physics_version", "initialization_method",
        "mip_era", "grid_label", "variant_label"
    ]
    n_members = length(model_names)
    members = Vector{String}(undef, n_members)
    for key in filter(k -> k in keys_model_ids, keys(meta))
        val = meta[key]
        if isa(val, String)
            meta_subdict[key] = repeat([val], outer=n_members)
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
            meta_subdict["physics_version"][indices_cmip5]
        )
        members[indices_cmip5] =  map(
            x -> join(x, MODEL_MEMBER_DELIM, MODEL_MEMBER_DELIM),
            zip(model_names[indices_cmip5], variants)
        )
        @debug "For CMIP5, full model names dont include grid."
    end
    if !isempty(indices_cmip6)
        variants = meta_subdict["variant_label"][indices_cmip6]
        grids = meta_subdict["grid_label"][indices_cmip6]
        members[indices_cmip6] = map(
            x->join(x, MODEL_MEMBER_DELIM, "_"),
            zip(model_names[indices_cmip6], variants, grids)
        );
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
    return result
end

"""
    subsetPaths(paths::Vector{String}, shared_models::Vector{String})

Return subset of paths in `paths` that point to models in `shared_models`.

# Arguments:
- `paths`: paths to model data
- `shared_models`: have form 'modelname#memberID[_grid]
"""
function subsetPaths(paths::Vector{String}, shared_models::Vector{String})
    # update metadata paths too, keep only those that contain a model 
    # in shared_models (NOTE: this doesnt work if the filename does 
    # not contain the model name, which it should though)
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
    subsetModelData(data::DimArray, shared_models::Vector{String})

Return a subset of DimArray `data` that contains only data from the models 
specified in `shared_models`. Takes care of metadata.

# Arguments:
- `data`:
- `shared_models`: 
"""
function subsetModelData(data::DimArray, shared_models::Vector{String})
    dim_symbol = hasdim(data, :member) ? :member : :model
    indices = findall(m -> m in shared_models, collect(dims(data, dim_symbol)))
    @assert length(indices) == length(shared_models)
    if dim_symbol == :model
        data = data[model = indices]
    else
        data = data[member = indices]
    end
    # also adjust the metadata
    attributes = filter(k -> data.metadata[k] isa Vector, keys(data.metadata))
    for key in attributes
        data.metadata[key] = data.metadata[key][indices]
    end
    return data
end


"""
    loadPreprocData(meta_data::MetaData, is_model_data::Bool)

# Arguments:
- `path_data`: base path to directory that contains subdirectory with name
of the respective experiment (see TODO for assumed data structure)
- `is_model_data`: observational and model data loaded seperately,
if true modelData, else observational data

# Returns: instance of `Data` or nothing
"""
function loadPreprocData(meta::MetaData, is_model_data::Bool=true)
    data = []
    meta_dict = Dict{String, Any}()
    n_files = length(meta.paths)
    source_names = repeat([""], outer = n_files)
    for (i, file) in enumerate(meta.paths)
        @debug "processing file.." * file
        parts = splitpath(file)
        filename = split(parts[end], ".nc")[end-1]
        climVar = split(parts[end-1], "_")[1]

        ds = NCDataset(file)
        dsVar = ds[climVar]
        # if climVar == "amoc"
        #     dsVar = ds["msftmz"]
        # end
        attributes = merge(Dict(deepcopy(dsVar.attrib)), Dict(deepcopy(ds.attrib)))
        if climVar == "msftmz"
            sector = get(ds, "sector", nothing)
            if !isnothing(sector)
                attributes = merge(attributes, Dict(deepcopy(sector.attrib)))
            end
        end
        # add mip_era for models since it is not provided in CMIP5-models
        name = filename
        if is_model_data
            model_key = getCMIPModelsKey(Dict(ds.attrib))
            get!(attributes, "mip_era", "CMIP5")
            name = ds.attrib[model_key] # just the model name, e.g. ACCESS1-0 (not the member's id)
            if occursin("_", name)
                @warn "model name as read from metadata of stored .nc file contains underscore: $name; from $file"
            end
        end
        source_names[i] = name
        # update metadata-dictionary for all processed files with the
        # metadata from the current file
        for key in keys(attributes)
            values = get!(meta_dict, key, repeat(Union{Missing, Any}[missing], outer=n_files))
            values[i] = attributes[key]
        end

        dimension_names = dimnames(dsVar)
        dimensions = []
        for d in dimension_names
            if d in ["bnds", "string21"]
                continue
            end
            if d == "time"
                times = map(x -> NCDatasets.DateTimeStandard(
                    Dates.year(x), Dates.month(x), Dates.day(x)
                    ),
                    dsVar[d][:]
                )
                push!(dimensions, Dim{Symbol(d)}(collect(times)))
            else
                push!(dimensions, Dim{Symbol(d)}(collect(dsVar[d][:])))
            end
        end
        push!(data, DimArray(Array(dsVar), Tuple(dimensions)))
    end
    if length(data) > 0
        #dimData = cat(data..., dims=3) # way too slow!
        # all of the preprocessed model data assumed to have the same grid!
        # Preallocate an n-dimensional array of respective size
        size_dims = size(data[1])
        n = length(data)
        raw_data = Array{eltype(Array(data[1]))}(undef, size_dims..., n)
        s = repeat([:], length(size_dims))
        names = collect(source_names)
        updateMetadata!(meta_dict, names, is_model_data)
        # sort the data
        sort_indices = sortperm(names)
        for idx in sort_indices
            raw_data[s..., idx] = Array(data[idx])
        end
        dimData = DimArray(
            raw_data,
            (dims(data[1])..., Dim{:source}(names[sort_indices]))
        )
        dimData = rebuild(dimData; metadata = meta_dict)
        if is_model_data
            indices = sortperm(dimData.metadata["member_names"])
            for idx in indices
                raw_data[s..., idx] = Array(data[idx])
            end
            dimData = DimArray(
                raw_data,
                (dims(data[1])..., Dim{:member}(dimData.metadata["member_names"][indices]))
            )
            dimData = rebuild(dimData; metadata = meta_dict) 

            # Sanity checks that no dataset exists more than once
            members = dims(dimData, :member)
            if length(members) != length(unique(members))
                duplicates = [m for m in members if sum(members .== m) > 1]
                @warn "Some datasets appear more than once" duplicates
            end
        end
        return Data(meta = meta, data = dimData)
    else
        return nothing
    end
end


"""
    loadDataFromMetadata(
        meta_data::Vector{MetaData},
        is_model_data::Bool,
        only_shared_models::Bool
    )

# Arguments:
- `meta_data`:
- `is_model_data`: true for model data, false for observational data.
- `only_shared_models`: if true, only data loaded for models present for all
ids, i.e. for all variable+statistic+experiment combinations
"""
function loadDataFromMetadata(
    meta_data::Dict{String, MetaData},
    is_model_data::Bool,
    only_shared_models::Bool
)
    results = Dict{String, Data}()
    for (id, meta) in meta_data
        results[id] = loadPreprocData(meta, is_model_data)
    end
    # The following shouldnt be necessary as metadata should already only contain 
    # shared models, but the following code retrieves shared models not from the 
    # filenames but from dimension names, which originate
    # from the metadata in the stored files. If there was a problem 
    # (e.g. FGOALS-g2 (correct, like filename) vs. FGOALS_g2 (unlike filename) for tos data)
    # this data will be excluded here. If it's not, dimension names might not be identical,
    # which seems to be problematic.
    # TODO: add checks/warnings here
    if is_model_data && only_shared_models
        shared_models = getSharedModels(results, :member)
        if isempty(shared_models)
            @warn "No models shared across all loaded data"
        end
        shared_results = Dict{String, Data}()# Vector{Data}()
        for (id, result) in results
            take_subset = sort(Array(dims(result.data, :member))) != sort(Array(shared_models))
            shared_data = take_subset ? subsetModelData(result.data, shared_models) : result.data
            shared_paths = take_subset ? subsetPaths(result.meta.paths, shared_models) : result.meta.paths
            meta = MetaData(id = result.meta.id, attrib = result.meta.attrib, paths = shared_paths)
            shared_results[id] = Data(meta = meta, data = shared_data)
        end
        results = shared_results
    end
    @debug "loaded data: $(map(x -> x.meta, values(results)))"
    @debug "filtered for shared models across all loaded data: " only_shared_models
    return results
end


"""
    getSharedModels(model_data::Vector{Data}, dim::Symbol)

Return only data of models shared across all datasets. Assumes that no id is 
shared across elements of `model_data`.

# Arguments
- `model_data`: model predictions for a set of different variable+experiment combinations
- `dim`: dimension name referring to level of model predictions; e.g., 
'member' or 'model'
"""
function getSharedModels(model_data::Dict{String, Data}, dim::Symbol)
    shared_models =  nothing 
    for (_, data) in model_data
        models = collect(dims(data.data, dim))
        if isnothing(shared_models)
            shared_models = models
        else
            shared_models = intersect(shared_models, models)
        end
    end
    return shared_models
end



"""
    getModelIDsFromPaths(all_paths::Vector{String})

For every path in `all_paths` returns a string of the form modelname#memberID[_grid]
that identifies the corresponding model (on the level of model members!).

# Arguments:
- `all_paths`:
"""
function getModelIDsFromPaths(all_paths::Vector{String})
    all_filenames = map(p -> split(basename(p), "_"), vcat(all_paths...))
    # model names are at predefined position in filenames (ERA_name_mip_exp_id_variable[_grid].nc)
    all_members = Vector{String}()
    for fn_parts in all_filenames
        model = join(fn_parts[[2, 5]], MODEL_MEMBER_DELIM)
        if fn_parts[1] != "CMIP5"
            model = join([model, fn_parts[7]], "_")
            model = split(model, ".nc")[1]
        end
        push!(all_members, model)
    end
    return unique(all_members)
end


"""
    searchModelInPaths(model::String, paths::Vector{String})

# Arguments:
- `model_id`: string with form: modelname#memberID[_grid]
- `paths`:
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
        meta_data::Dict{String, MetaData}, all_models::Vector{String}
    )
     
    Return a vector of models from `meta_data` that appear in `all_models`.

# Arguments:
- `meta_data`:
- `all_models`:
"""
function getSharedModelsFromPaths(
    meta_data::Dict{String, MetaData}, all_models::Vector{String}
)
    indices_shared = []
    for (idx, model) in enumerate(all_models)
        is_found = false
        for (_, meta) in meta_data
            is_found = searchModelInPaths(model, meta.paths)
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

Return the respective key to retrieve model names in CMIP6 and CMIP5 data.
"""
function getCMIPModelsKey(meta::Dict)
    attributes = keys(meta)
    if "source_id" in attributes
        return "source_id"
    elseif "model_id" in attributes
        return "model_id"
    else
        msg = "Only CMIP6/5 supported (model name keys: source_id/model_id)"
        throw(ArgumentError(msg))
    end
end


"""
    getUncertaintyRanges(data::DimArray, w::DimArray; quantiles=[0.167, 0.833]})

# Arguments:
- `data`: has dimensions 'time', 'model'
- `w`: has dimension 'model', sum up to 1
- `quantiles`: Vector with two entries (btw. 0 and 1) [lower_bound, upper_bound]
"""
function getUncertaintyRanges(data::DimArray, w::DimArray; quantiles=[0.167, 0.833])
    unweightedRanges = []
    weightedRanges = []
    for t in dims(data, :time)
        lower, upper = computeInterpolatedWeightedQuantiles(
            quantiles, Array(data[time = At(t)]); weights=w
        )
        push!(weightedRanges, [lower, upper])
        lower, upper = computeInterpolatedWeightedQuantiles(
            quantiles, Array(data[time = At(t)])
        )
        push!(unweightedRanges, [lower, upper])
    end

    return (weighted=weightedRanges, unweighted=unweightedRanges)
end


function showDataPaths(data::Union{Data, Vector{Data}})
    data = isa(data, Vector) ? data : [data]
    map(x->println(x.meta), data);
    return nothing
end