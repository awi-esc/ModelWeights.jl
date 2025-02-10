using DimensionalData
using NCDatasets
using YAML
using Setfield

"""
    getUniqueMemberIds(meta::Dict, model_names::Vector{String})
Combine name of the model with the id of the respective model members. The
returned vector has length of the total number of model members. It contains 
for every model the unique ids of its members.

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
    return vcat(result...)
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
- `data`: must have dimension 'member' or 'model'
- `shared_models`: vector of models, which can either be on level of models 
or members of models
"""
function subsetModelData(data::DimArray, shared_models::Vector{String})
    dim_symbol = hasdim(data, :member) ? :member : :model
    dim_names = collect(dims(data, dim_symbol))
    if dim_symbol == :member
        models = map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), shared_models)
        # if shared_models is on the level of models, the following should be empty
        # otherwise, nothing is filtered out, and members is the same as shared_models 
        members = filter(x -> !(x in models), shared_models)    
        if !isempty(members) # shared models on level of members
            indices = findall(m -> m in members, dim_names)
        else
            # important not to use dim_names here, since e.g. model=AWI would be found in dim_names where model is actually AWI-X for instance
            models_data = getModelsFromMemberIDs(dim_names) # NOTE: should yield same: models_data = data.metadata["model_names"]
            indices = findall(m -> m in models, models_data)
        end
        data = data[member = indices]
    else 
        indices = findall(m -> m in shared_models, dim_names)
        data = data[model = indices]
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
                #paths = [filter(x -> occursin(split(dup, MODEL_MEMBER_DELIM)[1], x), meta.paths) for dup in duplicates]
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
        is_model_data::Bool
    )

# Arguments:
- `meta_data`:
- `is_model_data`: true for model data, false for observational data.
"""
function loadDataFromMetadata(
    meta_data::Dict{String, MetaData},
    is_model_data::Bool
)
    results = Dict{String, Data}()
    for (id, meta) in meta_data
        # loads data at level of model members
        @info "load $id"
        results[id] = loadPreprocData(meta, is_model_data)
    end
    @debug "loaded data: $(map(x -> x.meta, values(results)))"
    @debug "filtered for shared models across all loaded data: " level
    return results
end


"""
    getSharedModelsFromDimensions(model_data::Vector{Data}, dim::Symbol)

Return only data of models shared across all datasets. Assumes that no id is 
shared across elements of `model_data`.

# Arguments
- `model_data`: model predictions for a set of different variable+experiment combinations
- `dim`: dimension name referring to level of model predictions; e.g., 
'member' or 'model'
"""
function getSharedModelsFromDimensions(model_data::Dict{String, Data}, dim::Symbol)
    shared_models =  nothing 
    for (_, data) in model_data
        if hasdim(data.data, :model) && dim == :member
            @warn "Model data has dimension :model, but subset of models requested to be taken wrt to level of model members! Level of models is used instead."
            models = collect(dims(data.data, :model))
        elseif hasdim(data.data, :member) && dim == :model
            members = collect(dims(data.data, :member))
            models = unique(getModelsFromMemberIDs(members))
        else
            models = collect(dims(data.data, dim))
        end
        if isnothing(shared_models)
            shared_models = models
        else
            shared_models = intersect(shared_models, models)
        end
    end
    return shared_models
end



"""
    getMemberIDsFromPaths(all_paths::Vector{String})

For every path in `all_paths` returns a string of the form modelname#memberID[_grid]
that identifies the corresponding model (on the level of model members!).

# Arguments:
- `all_paths`:
"""
function getMemberIDsFromPaths(all_paths::Vector{String})
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

Return the respective key to retrieve model names in CMIP6 ('source_id') and CMIP5 ('model_id') data.
If both keys are present, 'source_id' used in CMIP6 models is returned, if none is present, throw 
ArgumentError.

# Arguments:
- `meta`:
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


function reduceMetaDataSharedModels!(
    meta_data::Dict{String, MetaData}, level_shared_models::Union{LEVEL, Nothing}
)
    all_paths = map(p -> p.paths, values(meta_data))
    all_models = getMemberIDsFromPaths(vcat(all_paths...))
    if level_shared_models == MODEL
        all_models = unique(getModelsFromMemberIDs(all_models))
    end
    shared_models = getSharedModelsFromPaths(meta_data, all_models)
    if isempty(shared_models)
        @warn "No models shared across data!"
    end
    for (id, meta) in meta_data
        shared_paths = subsetPaths(meta.paths, shared_models)
        meta_new = @set meta.paths = shared_paths
        meta_data[id] = meta_new
    end
end

function getModelsFromMemberIDs(members::Vector{String})
    return map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), members)
end


function getPhysicsFromMembers(members::Vector{String})
    regex = r"(p\d+)(f\d+)?(_.*)?$"
    return String.(map(x -> match(regex, x).captures[1], members))
end


"""
    alignPhysics(data::Dict{String, Data}, members::Vector{String})

Return new data dictionary with only the models retained that share the same physics 
as the respective model's members in `members` have.

# Arguments:
- `data`: 
- `members`:
"""
function alignPhysics(
    data::Dict{String, Data}, members::Vector{String};
    level_shared_models::Union{LEVEL, Nothing} = nothing
)
    models = unique(getModelsFromMemberIDs(members))
    for model in models
        # retrieve allowed physics as in members 
        member_ids = filter(m -> startswith(m, model * MODEL_MEMBER_DELIM), members)
        physics = getPhysicsFromMembers(member_ids)
        for (id, model_data) in data
            ds = model_data.data
            # filter data s.t. of current model only members with retrieved physics are kept
            model_indices = findall(x -> startswith(x, model * MODEL_MEMBER_DELIM), Array(dims(ds, :member)))
            indices_out = filter(x -> !(ds.metadata["physics"][x] in physics), model_indices)
            if !isempty(indices_out)
                indices_keep = filter(x -> !(x in indices_out), 1:length(dims(ds, :member)))
                members_kept = ds.metadata["member_names"][indices_keep]
                meta = model_data.meta
                meta_updated = @set meta.paths = subsetPaths(meta.paths, members_kept)
            
                data[id] = Data(meta = meta_updated, data = subsetModelData(ds, members_kept))
            end
        end
    end
    # after having subset the data wrt the physics, make sure that level of 
    # shared models as it was before if wanted
    if !isnothing(level_shared_models)
        meta_data = Dict{String, MetaData}()
        for k in collect(keys(data))
            meta_data[k] = data[k].meta 
        end
        reduceMetaDataSharedModels!(meta_data, level_shared_models)
        all_paths = reduce(vcat, map(k -> meta_data[k].paths, collect(keys(data))))
        all_members = getMemberIDsFromPaths(all_paths)
        for (id, model_data) in data
            data[id] = Data(
                meta = meta_data[id], 
                data = subsetModelData(model_data.data, all_members)
            )
        end
    end
    return data
end

""" 
    summarizeEnsembleMembersVector(
        data::DimArray, updateMeta::Bool; fn::Function=Statistics.mean
)

For each model and variable (if several given), compute a summary statistic 
(default: mean) across all members of that model. Instead of 'member', the 
returned DimArray has dimension 'model'.

# Arguments:
- `data`: a DimArray with at least dimension 'model'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
Set to false if vectors refer to different variables for instance. 
"""
function summarizeEnsembleMembersVector(
    data::DimArray, updateMeta::Bool; fn::Function=Statistics.mean
)
    data = setLookupsFromMemberToModel(data, ["member"])
    grouped = groupby(data, :model=>identity);
    models = String.(collect(dims(grouped, :model)))
    averages = map(entry -> mapslices(x -> fn(skipmissing(x)), entry, dims=:model), grouped)
    combined = cat(averages..., dims=(Dim{:model}(models)));
    combined = replace(combined, NaN => missing)

    meta = updateMeta ? updateGroupedDataMetadata(data.metadata, grouped) : data.metadata
    combined = rebuild(combined; metadata = meta);
    l = Lookups.Categorical(
        sort(models);
        order=Lookups.ForwardOrdered()
    )
    combined = combined[model=At(sort(models))]
    combined = DimensionalData.Lookups.set(combined, model=l)

    return combined
end


"""
    averageEnsembleMembers!(data::Dict{String, Data})

Take average for all members of each model.

# Arguments:
- `data`: 
"""
function averageEnsembleMembers!(data::Dict{String, Data})
    for k in keys(data)
        updated_arr = summarizeEnsembleMembersVector(data[k].data, true)
        current_data = data[k]
        data[k] = @set current_data.data = updated_arr 
    end
    return nothing
end


"""
    getGlobalMeans(data::DimArray)

Compute area-weighted globalMeans across longitudes and latitudes for each 
model. Missing data is accounted for in the area-weights. 

# Arguments:
- `data`: DimArray with at least dimensions lon, lat and possibly member or model.

# Return: a DimArray of size 'number models in data' x 1 containing area-weighted
global means for each model. 
"""
function getGlobalMeans(data::DimArray)
    longitudes = Array(dims(data, :lon))
    latitudes = Array(dims(data, :lat))
    masks = ismissing.(data)

    dimension = hasdim(data, :member) ? :member : hasdim(data, :model) ? :model : nothing
    if !isnothing(dimension)
        models = Array(dims(data, dimension))
        global_means = DimArray(zeros(length(models)), (Dim{dimension}(models)))
        for model in models
            if dimension == :model
                mask = masks[model = At(model)]
                area_weights = computeAreaWeights(longitudes, latitudes; mask)
                global_means[model = At(model)] = Statistics.sum(data[model = At(model)] .* area_weights)
            else
                mask = masks[member = At(model)]
                area_weights = computeAreaWeights(longitudes, latitudes; mask)
                global_means[member = At(model)] = Statistics.sum(skipmissing(data[member = At(model)] .* area_weights))
            end
        end
    else 
        area_weights = computeAreaWeights(longitudes, latitudes; mask=masks)
        global_means = Statistics.sum(skipmissing(data .* area_weights))
    end
    return global_means
end




function compute_anomalies(orig::Data, ref::Data; save::Bool=true)
    dimension = hasdim(orig.data, :member) ? :member : :model
    if dims(orig.data, dimension) != dims(ref.data, dimension)
        throw(ArgumentError("Original and reference data must contain exactly the same models!"))
    end
    anomalies = orig.data .- ref.data
    # getPaths(data::Data) = map(x -> (base=joinpath(splitpath(x)[1:end-3]), 
    #                                 alias=splitpath(x)[end-2],
    #                                 diagnostic=splitpath(x)[end-1],
    #                                 fn = splitpath(x)[end]), 
    #                             data.meta.paths)
    # # if save 
    #     paths_refs = getPaths(ref)
    #     paths_data = getPaths(orig)
    # end
    return anomalies
end