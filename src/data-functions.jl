using DimensionalData

"""
    getUniqueModelIds(meta::Dict, key_model_names::String)

# Arguments:
- `meta`:
- `key_model_names`:
"""
function getUniqueModelIds(
    meta::Dict, #{String, Union{String, Array, Dict}},
    key_model_names::String
)
    mip_era = get(meta, "mip_era", get(meta, "project_id", ""))
    model_names = meta[key_model_names];
    if mip_era == "CMIP5"
        variants = buildCMIP5EnsembleMember(
            meta["realization"], 
            meta["initialization_method"], 
            meta["physics_version"],
            length(model_names)
        );
        models = map(
            x->join(x, MODEL_MEMBER_DELIM, MODEL_MEMBER_DELIM), 
            zip(model_names, variants)
        );
        @warn "For CMIP5, full model names dont include grid."
    else
        variants = meta["variant_label"];
        grids = meta["grid_label"]; 
        models = map(
            x->join(x, MODEL_MEMBER_DELIM, MODEL_MEMBER_DELIM), 
            zip(model_names, variants, grids)
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
    loadPreprocData(path_data_dir::String, included::Vector{String}=[])

    This is the new version!

# Arguments:
- `path_data_dir`: base path to directory that contains subdirectory with name
of the respective experiment (see TODO for assumed data structure)
- `included`: either observational or model data can be loaded at once, e.g. 
set to ["ERA5"] or ["CMIP5"] for ERA5-observational and CMIP5 model data respectively
"""
function loadPreprocData(path_data_dir::String, included::Vector{String}=[])
    if !isdir(path_data_dir)
        throw(ArgumentError(path_data_dir * " does not exist!"))
    end
    data = [];
    sources = [];
    meta = Dict{String, Union{String, Array, Dict}}();
    ncFiles = filter(
        x -> isfile(x) && endswith(x, ".nc"), 
        readdir(path_data_dir; join=true)
    );
    nbIgnored = 0;
    for (i, file) in enumerate(ncFiles)
        addFile = true;
        # only include files that contain all names given in 'included'
        if length(included) != 0
            if !all([occursin(name, file) for name in included])
                addFile = false;
            end
        end
        #println("processing file.." * file)
        if addFile
            parts = splitpath(file)
            filename = split(parts[end], ".nc")[end-1]
            climVar = split(parts[end-1], "_")[1]
            
            ds = NCDataset(file);
            dsVar = ds[climVar];
            # if climVar == "amoc"
            #     dsVar = ds["msftmz"]
            # else 
            # end
            
            attributes = merge(Dict(deepcopy(dsVar.attrib)), Dict(deepcopy(ds.attrib)));
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

            # update metadata
            for key in keys(attributes)
                values = get!(meta, key, []);
                # fill up vector
                n = i - nbIgnored - length(values) - 1;
                for _ in range(1, n)
                    push!(values, missing)
                end
                push!(values, attributes[key]);
            end

            name = get(ds.attrib, "source_id", get(ds.attrib, "model_id", ""));
            if !isempty(name)
                push!(sources, name);
            else # observational data does not have model name
                push!(sources, filename)
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
        dimData = DimArray(raw_data, (DimensionalData.dims(data[1])..., Dim{:model}(sources)))
        updateMetadata!(meta, dimData);
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
        subset::Dict=Dict(),
    )

Loads the data from the config files located at 'paths_to_config_dir'. For necessary
structure of config files, see TODO. For each variable, experiment, statistic
and timerange (task) a different DimArray is loaded.

# Arguments:
- `base_path`:  if dir_per_var is true, path to directory that contains one or
more subdirectories that each contains a directory 'preproc' with the
preprocessed data. If dir_per_var is false, base_path is the path to a directory
that contains the 'preproc' subdirectory.
- `config_path`: path to directory that contains one or more yaml config 
files with the following structure: TODO
- `dir_per_var`: if true, one subdirectory for data of each climate variable
- `subset`: dictionary specifying the subset of data to be loaded, has keys
'variables', 'statistics', 'aliases', each mapping to a vector of Strings
"""
function loadData(
    base_path::String,
    config_path::String;
    dir_per_var::Bool=true,
    common_models_across_vars=false,
    subset::Dict{String, Vector{String}}=Dict{String, Vector{String}}()
)
    ids = buildDataIDsFromConfigs(config_path)    
    applyDataConstraints!(ids, subset)

    model_data = Dict{String, DimArray}()
    obs_data = Dict{String, DimArray}()
    for id in ids
        # if dir_per_var is true, directory at base_path has subdirectories, 
        # one for each variable (they must end with _ and the name of the variable),
        # otherwise config_path is the path to the directory that contains 
        # a subdirectory 'preproc'
        path_to_subdirs = [base_path]
        if dir_per_var
            path_to_subdirs = filter(isdir, readdir(base_path, join=true))
            filter!(x -> endswith(x, "_" * id.variable), path_to_subdirs)
            if length(path_to_subdirs) > 1
                @warn "There are several subdirectories for variable " * var * " and experiment " * exp
                path_to_subdirs = path_to_subdirs[1:1]
                @warn "First is chosen: " path_to_subdirs[1]
            end
        end
        for path_dir in path_to_subdirs
            path_data_dir = joinpath(
                path_dir, "preproc", id.task, join([id.variable, id.statistic], "_")
            )
            if !isdir(path_data_dir)
                continue
            end
            # the name of the time period is needed for loading the data, but 
            # as it is arbitrary (e.g. historical1 ~ 1950-1981), we do not 
            # included in the data ids
            #base_id = join([id.statistic, id.variable, id.timerange], "_")
            base_id = join(split(id.key, "_")[2:end], "_")
            #id_wit_exp =  join([id.exp, base_id], "_")
            data = loadPreprocData(path_data_dir, ["CMIP"])
            if !isnothing(data)
                #data = convert(DimArray, data)
                model_data[id.key] = data
            end
            # TODO: don't hard code name of observational dataset!
            name_obs_data = "ERA5"
            obs = loadPreprocData(path_data_dir, [name_obs_data])
            if !isnothing(obs)
                #obs = convert(Dict{String, DimArray}, obs)
                obs_data[join([name_obs_data, base_id], "_")] = obs
            end
        end
    end
    @info "The following data was found and loaded: " keys(model_data)

    if common_models_across_vars
        model_data = getCommonModelsAcrossVars(model_data, ids)
    end
    return Data(
        base_path = base_path,
        ids = ids,
        models = model_data,
        obs = obs_data
    )
end


"""
    getCommonModelsAcrossVars(modelData::Dict{String, DimArray}, ids::Vector{DataID})

Return only those models for which there is data for all variables.

# Arguments
- `modelData`:
- `ids`:
"""
function getCommonModelsAcrossVars(modelData::Dict{String, DimArray}, ids::Vector{DataID})
    data_all = deepcopy(modelData);
    variables = map(id -> id.variable, ids)
    shared_models =  nothing
    for var in variables
        modelDict = filter(((k,v),)-> occursin(var, k), modelData)
        if length(modelDict) > 1
            @warn "more than one dataset for computing getCommonModelsAcrossVars, first is taken!"
        end
        data_var = first(values(modelDict))
        meta = data_var.metadata;
        models = getUniqueModelIds(meta, getCMIPModelsKey(meta));
        data_all[first(keys(modelDict))].metadata["full_model_names"] = models
        if isnothing(shared_models)
            shared_models = models
        else
            shared_models = intersect(shared_models, models)
        end
    end

    for id in map(id -> id.key, ids)
        data_all[id] = getModelSubset(data_all[id], shared_models);
    end

    return data_all
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
        msg = "Only CMIP6/CMIP5 supported with model name keys: source_id/model_id!"
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