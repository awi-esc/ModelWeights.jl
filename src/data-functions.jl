using DimensionalData

"""
    getUniqueModelIds(meta::Dict, model_names::Vector{String})


# Arguments:
- `meta`: for CMIP5 data, must have keys: 'mip_era', 'realization', 'initialization_method',
'physics_version'. For CMIP6 data must have keys: 'variant_label', 'grid_label'
- `model_names`:
"""
function getUniqueModelIds(
    meta::Dict,
    model_names::Vector{String}
)
    meta_subdict = Dict{String, Vector}()
    keys_model_ids = [
        "realization", "physics_version", "initialization_method", 
        "mip_era", "grid_label", "variant_label"
    ]
    n = length(model_names)
    for key in filter(k -> k in keys_model_ids, keys(meta))
        val = meta[key]
        if isa(val, String)
            meta_subdict[key] = repeat([val], outer=n)
        else
            meta_subdict[key] = val
        end
    end
    mip_eras = meta_subdict["mip_era"]
    indices_cmip5 = findall(x -> !ismissing(x) && x == "CMIP5", mip_eras)
    indices_cmip6 = findall(x -> !ismissing(x) && x == "CMIP6", mip_eras)
    models = Vector{String}(undef, length(model_names))

    if !isempty(indices_cmip5)
        variants = buildCMIP5EnsembleMember(
            meta_subdict["realization"][indices_cmip5], 
            meta_subdict["initialization_method"][indices_cmip5], 
            meta_subdict["physics_version"][indices_cmip5]
        )
        models[indices_cmip5] =  map(
            x -> join(x, MODEL_MEMBER_DELIM, MODEL_MEMBER_DELIM), 
            zip(model_names[indices_cmip5], variants)
        )
        @debug "For CMIP5, full model names dont include grid."
    end
    if !isempty(indices_cmip6)
        variants = meta_subdict["variant_label"][indices_cmip6]
        grids = meta_subdict["grid_label"][indices_cmip6]
        models[indices_cmip6] = map(
            x->join(x, MODEL_MEMBER_DELIM, "_"), 
            zip(model_names[indices_cmip6], variants, grids)
        );
    end
    return models
end

"""
    getIndicesMapping(names)

# Arguments:
- `names`:
"""
function getIndicesMapping(names)
    mapping = Dict();
    for name in names
        mapping[name] = findall(x -> x == name, names)
    end
    return mapping
end

"""
    getModelSubset(data::Dict{String, DimArray}, shared_models::Vector{String})

Return data in 'data' only from the models specified in `shared_models`. 
Takes care of metadata.
"""
function getModelSubset(data::DimArray, shared_models::Vector{String})
    indices = findall(m -> m in shared_models, data.metadata["full_model_names"]);
    @assert length(indices) == length(shared_models)
    keepMetadataSubset!(data.metadata, indices);
    data = data[model = indices];
    data.metadata["indices_map"] = getIndicesMapping(data.metadata["ensemble_names"])
    return data
end


"""
    getSharedMetadataAndModelNames(metadata::Dict)

Return a new metadata dictionary which contains all attributes that were
identical across models (therefore these are single values, not Vectors). 
Further, model names and full_model_names are added. When combining this new 
metadata dict with another, e.g. when combining data for different variables, 
these must be identical (which is checked in function appendValuesDicts).
Ignores indices_map (mapping models to indices in metadata arrays), since the 
shared data may be combined differently.
"""
function getSharedMetadataAndModelNames(metadata::Dict)
    models_key = getCMIPModelsKey(metadata);
    meta_shared = filter(((k,v),) -> isa(v, String), metadata);
    meta_shared["full_model_names"] = deepcopy(metadata)["full_model_names"];
    meta_shared[models_key] = unique(metadata[models_key]);
    return meta_shared
end


"""
    loadPreprocData(
        path_data_dir::String;
        included_all::Vector{String}=Vector{String}(),
        included_any::Vector{String}=Vector{String}(),
        isModelData::Bool=true
    )

# Arguments:
- `path_data_dir`: base path to directory that contains subdirectory with name
of the respective experiment (see TODO for assumed data structure)
- `included_all`: only data is loaded whose filenames contain ALL elements within 'included_all'
- `included_any`: only data is loaded whose filenames containy ANY element within 'included_any'
- `isModelData`: observational and model data loaded seperately, if true modelData, else observational data

# Returns: DimArray or nothing
"""
function loadPreprocData(
    path_data_dir::String; 
    included_all::Vector{String}=Vector{String}(),
    included_any::Vector{String}=Vector{String}(),
    isModelData::Bool=true
)
    if !isdir(path_data_dir)
        throw(ArgumentError(path_data_dir * " does not exist!"))
    end
    # set default values for model data
    if isempty(included_all)
        if isModelData
            included_all = ["CMIP"]
        else
            included_all = ["ERA5"]
        end
    end
    data = []
    meta = Dict{String, Any}()# Dict{String, Union{String, Vector, Dict}}()
    ncFiles = filter(
        x -> isfile(x) && endswith(x, ".nc"), 
        readdir(path_data_dir; join=true)
    )
    nbIgnored = 0
    n_files = length(ncFiles)
    source_names = repeat(Union{Missing, String}[missing], outer = n_files)
    for (i, file) in enumerate(ncFiles)
        addFile = true;
        # constrain files that will be loaded
        if length(included_all) != 0
            if !all([occursin(name, file) for name in included_all])
                addFile = false;
            end
        end
        if length(included_any) != 0
            if !any([occursin(name, file) for name in included_any])
                addFile = false;
            end
        end
        if addFile
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
                    attributes = merge(attributes, Dict(deepcopy(sector.attrib)));
                end
            end
            if warnIfFlawedMetadata(attributes, filename)
                nbIgnored += 1;
                continue
            end
            # add mip_era for models since it is not provided in CMIP5-models
            name = ""
            if isModelData
                model_key = getCMIPModelsKey(Dict(ds.attrib))
                get!(attributes, "mip_era", "CMIP5")
                name = ds.attrib[model_key]
            else 
                name = filename
            end
            source_names[i] = name        
            # update metadata-dictionary for all processed files with the 
            # metadata from the current file
            for key in keys(attributes)
                values = get!(meta, key, repeat(Union{Missing, Any}[missing], outer=n_files));
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
                    );
                    push!(dimensions, Dim{Symbol(d)}(collect(times)));
                else
                    push!(dimensions, Dim{Symbol(d)}(collect(dsVar[d][:])));
                end
            end      
            push!(data, DimArray(Array(dsVar), Tuple(dimensions)));
        end
    end
    if length(data) > 0
        #dimData = cat(data..., dims=3) # way too slow!
        # all of the preprocessed model data assumed to have the same lon,lat grid
        # Preallocate a 3D array of respective size
        size_dims = size(data[1])
        n = length(data)
        raw_data = Array{eltype(parent(data[1]))}(undef, size_dims..., n)
        if length(size_dims) == 1
            for idx in 1:n
                raw_data[:,idx] = parent(data[idx])
            end
        elseif length(size_dims) == 2
            for idx in 1:n
                raw_data[:, :, idx] = parent(data[idx])  # Extract and assign the underlying data
            end
        else 
            throw(ArgumentError("More than 3 data dimensions not supported."))
        end
        updateMetadata!(meta, source_names, isModelData);
        dimData = DimArray(raw_data, (dims(data[1])..., Dim{:model}(collect(skipmissing(source_names)))))
        dimData = rebuild(dimData; metadata = meta);

        # set model names in model dimension to full model name which was 
        # built in updateMetadata! (cannot be done before since missing 
        # files etc. have to be accounted for first)
        if !isempty(get(meta, "full_model_names", []))
            dimData = set(dimData, :model => meta["full_model_names"])
        end
        return dimData
    else
        return nothing
    end
end



"""
    loadData(
        base_path::String, 
        config_path::String;
        dir_per_var::Bool=true,
        isModelData::Bool=true,
        common_models_across_vars::Bool=false,
        subset::Dict=Dict(),
    )

Loads the data from the config files located at 'config_path'. For necessary
structure of config files, see TODO. For each variable, experiment, statistic
and timerange (alias) a different DimArray is loaded.

# Arguments:
- `base_path`:  if dir_per_var is true, path to directory that contains one or
more subdirectories that each contains a directory 'preproc' with the
preprocessed data. If dir_per_var is false, base_path is the path to a directory
that contains the 'preproc' subdirectory.
- `config_path`: path to directory that contains one or more yaml config 
files with the following structure: TODO
- `dir_per_var`: if true, directory at base_path has subdirectories, one for
each variable (they must end with _ and the name of the variable), otherwise
base_path is the path to the directory that contains a subdirectory 'preproc'
- `isModelData`: set true for CMIP5/6 data, false for observational data
- `common_models_across_vars`:
- `subset`: dictionary specifying the subset of data to be loaded, has keys
'variables', 'statistics', 'aliases', each mapping to a vector of Strings
"""
function loadData(
    base_path::String,
    config_path::String;
    dir_per_var::Bool=true,
    isModelData::Bool=true,
    common_models_across_vars::Bool=false,
    subset::Dict{String, Vector{String}}=Dict{String, Vector{String}}()
)
    ids = buildDataIDsFromConfigs(config_path)
    applyDataConstraints!(ids, subset)

    data_all = Dict{String, DimArray}()
    for id in ids
        # if dir_per_var is true, directory at base_path has subdirectories, 
        # one for each variable (they must end with _ and the name of the variable),
        # otherwise base_path is the path to the directory that contains 
        # a subdirectory 'preproc'
        path_to_subdirs = [base_path]
        if dir_per_var
            path_to_subdirs = filter(isdir, readdir(base_path, join=true))
            filter!(x -> endswith(x, "_" * id.variable), path_to_subdirs)
            if length(path_to_subdirs) > 1
                @warn "There are several subdirectories for given variable and experiment:" path_to_subdirs
                #path_to_subdirs = path_to_subdirs[1:1]
                #@warn "First is chosen: " path_to_subdirs[1]
            end
        end
        for path_dir in path_to_subdirs
            path_data_dir = joinpath(
                path_dir, "preproc", id.alias, join([id.variable, id.statistic], "_")
            )
            if !isdir(path_data_dir)
                @warn "$path_data_dir does not exist"
                continue
            end
            #print("processing...: " * path_data_dir)
            data = loadPreprocData(
                path_data_dir; 
                included_all = get(subset, "data_type", Vector{String}()),
                included_any = get(subset, "models", Vector{String}()),
                isModelData = isModelData
            )
            if !isnothing(data)
                previously_added_data = get(data_all, id.key, nothing)
                if isnothing(previously_added_data)
                    data_all[id.key] = data
                else
                    prev_models = collect(dims(previously_added_data, :model))
                    new_models = collect(dims(data, :model))
                    joint_data = cat(previously_added_data, data, dims=Dim{:model}(vcat(prev_models, new_models)))
                    joint_meta = joinMetadata(previously_added_data.metadata, data.metadata)
                    data_all[id.key] = rebuild(joint_data; metadata = joint_meta)
                end
            end
        end
    end
    result =  Data(base_path = base_path, ids = ids, data = data_all)
    @info "The following data was found and loaded: " result.ids
    if common_models_across_vars
        @info "only retain models shared across all variables"
        result = getCommonModelsAcrossVars(result)
    else
        alignIDsFilteredData!(result)
    end
    return result
end
    

"""
    getCommonModelsAcrossVars(modelData::Data)

Return only those models for which there is data for all variables.
TODO: add checks that model dimensions are identical everywhere!

# Arguments
- `modelData`
"""
function getCommonModelsAcrossVars(modelData::Data)
    data_all = deepcopy(modelData.data);
    variables = map(id -> id.variable, modelData.ids)
    shared_models =  nothing
    for var in variables
        modelDict = filter(((k,v),)-> occursin(var, k), modelData.data)
        # iterate over all combinations (of diagnostics/statistics) with current variable
        for (k, data_var) in modelDict
            meta = data_var.metadata;
            models = getUniqueModelIds(meta, Array(dims(data_var, :model)))
            data_all[k].metadata["full_model_names"] = models
            if isnothing(shared_models)
                shared_models = models
            else
                shared_models = intersect(shared_models, models)
            end
        end
    end

    for id in map(id -> id.key, modelData.ids)
        data_all[id] = getModelSubset(data_all[id], shared_models);
    end
    result = Data(base_path = modelData.base_path, ids = modelData.ids, data = data_all)
    alignIDsFilteredData!(result)
    return result
end


"""
    getCMIPModelsKey(meta::Dict)

Return the respective key to retrieve model names in CMIP6 and CMIP5 data.
"""
function getCMIPModelsKey(meta::Dict)
    attributes = keys(meta);
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
    saveWeights(
        weights::DimVector,
        target_dir::String;
        target_fn::String="weights.nc"
    )

# Arguments:
- `weights`:
- `target_dir`:
- `target_fn`:
"""
function saveWeights(
    weights::DimVector,
    target_dir::String;
    target_fn::String="weights.nc"
)
    if !isdir(target_dir)
        mkpath(target_dir)
    end
    path_to_target = joinpath(target_dir, target_fn);
    if isfile(path_to_target)
        msg1 = "File: " * path_to_target * " already exists!";
        path_to_target = joinpath(target_dir, join([getCurrentTime(), target_fn], "_"))
        msg2 = "Weights saved as: " * path_to_target
        @warn msg1 * msg2
    end
    ds = NCDataset(path_to_target, "c")

    models = dims(weights, :model)
    defDim(ds, "model", length(models))
    # Add a new variable to store the model names
    v_model_names = defVar(ds, "model", String, ("model",))
    v_model_names[:] = Array(models)

    # global attributes
    for (k, v) in weights.metadata
        if isa(v, Dict)
            ds.attrib[k] = [dk * "_" * string(dv) for (dk, dv) in v]
        elseif isa(v, String)
            ds.attrib[k] = v
        else
            ds.attrib[k] = Vector{String}(v)
        end
    end
    v = defVar(ds, "weight", Float64, ("model",))
    v[:] = Array(weights)
    close(ds)

    @info "saved data to " path_to_target
end


"""
    loadWeightsAsDimArray(path_to_file::String)
"""
function loadWeightsAsDimArray(path_to_file::String)
    data = NCDataset(path_to_file)
    models = Array(data["model"])
    arr = DimArray(
        Array(data["weight"]), 
        (Dim{:model}(models)), metadata = Dict(data.attrib)
    )
    return arr
end


"""
    getUncertaintyRanges(data::DimArray, w::DimArray; quantiles=[0.167, 0.833]})
    
# Arguments:
- `data`: has dimensions 'time', 'model'
- `w`: has dimension 'model', sum up to 1
- `quantiles`: Vector with two entries (btw. 0 and 1) [lower_bound, upper_bound]
"""
function getUncertaintyRanges(
    data::DimArray,
    w::DimArray;
    quantiles=[0.167, 0.833]
)
    unweightedRanges = [];
    weightedRanges = [];
    for t in dims(data, :time)
        lower, upper = computeInterpolatedWeightedQuantiles(
            quantiles, Array(data[time = At(t)]); weights=w
        );
        push!(weightedRanges, [lower, upper]);
        lower, upper = computeInterpolatedWeightedQuantiles(
            quantiles, Array(data[time = At(t)])
        );
        push!(unweightedRanges, [lower, upper]);
    end

    return (weighted=weightedRanges, unweighted=unweightedRanges)
end