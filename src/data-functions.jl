using DimensionalData
using NCDatasets

"""
    getUniqueModelIds(meta::Dict, model_names::Vector{String})

Combine name of the ensemble with the id of the respective ensemble members.

# Arguments:
- `meta`: for CMIP5 data, must have keys: 'mip_era', 'realization', 'initialization_method',
'physics_version'. For CMIP6 data must have keys: 'variant_label', 'grid_label'
- `model_names`: Vector of strings containing the names of the model ensembles.
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
            meta_subdict[key] = deepcopy(val)
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
    getIndicesMapping(names::Vector{String})

# Arguments:
- `names`: names of model members
"""
function getIndicesMapping(names::Vector{String})
    mapping = Dict();
    for name in names
        mapping[name] = findall(x -> x == name, names)
    end
    return mapping
end

"""
    subsetModelData(data::Dict{String, DimArray}, shared_models::Vector{String})

Return data in 'data' only from the models specified in `shared_models`. 
Takes care of metadata.
"""
function subsetModelData(data::DimArray, shared_models::Vector{String})
    dim_symbol = !hasdim(data, :model) ? :member : :model
    indices = findall(m -> m in shared_models, collect(dims(data, dim_symbol)))
    @assert length(indices) == length(shared_models)
    keepMetadataSubset!(data.metadata, indices);
    if dim_symbol == :model
        data = data[model = indices];
    else
        data = data[member = indices];
    end
    return data
end


"""
    getMetadataSharedAcrossModelsAndModelNames(metadata::Dict)

Return a new metadata dictionary which contains all attributes that were
identical across models (therefore these are single values, not Vectors). 
Further, 'model_names', 'member_names' and 'mip_era' are retained. When
combining this new metadata dict with another, e.g. when
combining data for different variables, these must be identical (which is 
checked in function appendValuesDicts).
"""
function getMetadataSharedAcrossModelsAndModelNames(metadata::Dict)
    meta_shared = filter(((k,v),) -> isa(v, String), metadata);    
    meta_shared["member_names"] = deepcopy(metadata)["member_names"];
    meta_shared["model_names"] = deepcopy(metadata)["model_names"]
    meta_shared["mip_era"] = deepcopy(metadata)["mip_era"]
    return meta_shared
end


"""
    loadPreprocData(
        path_data::String;
        subset::Dict{String, Vector{String}}=Dict{String, Vector{String}}(),
        is_model_data::Bool=true
    )

# Arguments:
- `path_data`: base path to directory that contains subdirectory with name
of the respective experiment (see TODO for assumed data structure)
- `subset`: dictionary with keys 'projects' and 'models'. Only data is loaded
whose filenames contain ANY of the strings mapped to. If not specified, default 
value for projects is ['CMIP'] when is_model_data is true, else ['ERA5']. 
- `is_model_data`: observational and model data loaded seperately, 
if true modelData, else observational data

# Returns: DimArray or nothing
"""
function loadPreprocData(
    path_data::String; 
    subset::Dict{String, Vector{String}}=Dict{String, Vector{String}}(),
    is_model_data::Bool=true
)
    if !isdir(path_data)
        throw(ArgumentError(path_data * " does not exist!"))
    end
    if isnothing(get(subset, "projects", nothing))
        subset["projects"] = is_model_data ? ["CMIP"] : ["ERA5"]
    end
    data = []
    meta = Dict{String, Any}()
    ncFiles = filter(
        x -> isfile(x) && endswith(x, ".nc"), 
        readdir(path_data; join=true)
    )
    nbIgnored = 0
    n_files = length(ncFiles)
    source_names = repeat(Union{Missing, String}[missing], outer = n_files)
    for (i, file) in enumerate(ncFiles)
        addFile = true;
        # constrain files that will be loaded
        if !any([occursin(name, file) for name in subset["projects"]])
            @debug "exclude $file because of projects subset"
            addFile = false;
        end
        model_constraints = get(subset, "models", Vector{String}())
        if !isempty(model_constraints)
            # model constraints may contain individual members
            model_member_constraints = map(x -> split(x, MODEL_MEMBER_DELIM), model_constraints)
            any_fullfilled = false
            for constraints in model_member_constraints
                constraint_ok = all([occursin(name, file) for name in constraints])
                if constraint_ok
                    any_fullfilled = true
                    break
                end
            end
            if !any_fullfilled
                @debug "exclude $file because of models subset"
                addFile = false
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
            # if warnIfFlawedMetadata(attributes, filename)
            #     nbIgnored += 1;
            #     continue
            # end
            # add mip_era for models since it is not provided in CMIP5-models
            name = ""
            if is_model_data
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
        dimData = DimArray(
            raw_data, 
            (dims(data[1])..., Dim{:source}(collect(skipmissing(source_names))))
        )
        updateMetadata!(meta, source_names, is_model_data);
        dimData = rebuild(dimData; metadata = meta);
        if is_model_data
            # set dimension names, member refers to unique model members, 
            # model refers to 'big combined model', part of member name,
            # but additionally saved in metadata["model_names"]
            dimData = set(dimData, :source => :member)
            dimData = set(dimData, :member => dimData.metadata["member_names"])
            # Sanity checks that no dataset exists more than once
            members = dims(dimData, :member)
            if length(members) != length(unique(members))
                duplicates = [m for m in members if sum(members .== m) > 1]
                @warn "Some datasets appear more than once" duplicates
            end
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
        is_model_data::Bool=true,
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
- `is_model_data`: set true for CMIP5/6 data, false for observational data
- `common_models_across_vars`:
- `subset`: dictionary specifying the subset of data to be loaded, has keys
'variables', 'statistics', 'aliases', each mapping to a vector of Strings
"""
function loadData(
    base_path::String,
    config_path::String;
    dir_per_var::Bool=true,
    is_model_data::Bool=true,
    common_models_across_vars::Bool=false,
    subset::Dict{String, Vector{String}}=Dict{String, Vector{String}}()
)
    ids = buildDataIDsFromConfigs(config_path)
    applyDataConstraints!(ids, subset)
    # further constraints wrt models and projects applied when loading data
    constraints = filter(((k,v),) -> k in ["models", "projects"] , subset)

    data_all = Dict{String, DimArray}()
    ids_all = Vector{DataID}()
    for id in ids
        path_to_subdirs = [base_path]
        # if dir_per_var is true, directory at base_path has subdirectories, 
        # one for each variable (they must contain '_VAR', e.g. '_tas'),
        # otherwise base_path is the path to the directory that contains 
        # a subdirectory 'preproc'
        if dir_per_var
            path_to_subdirs = filter(isdir, readdir(base_path, join=true))
            filter!(x -> occursin("_" * id.variable, x), path_to_subdirs)
            subdirs = get(subset, "subdirs", nothing)
            if !isnothing(subdirs)
                filter!(p -> any([occursin(name, p) for name in subdirs]), path_to_subdirs)
            end
            if length(path_to_subdirs) > 1
                fnames = map(basename, path_to_subdirs)
                @info "Data for variable $(id.key) considered from files:" fnames
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
                path_data_dir; subset = constraints, is_model_data = is_model_data
            )
            if !isnothing(data)
                previously_added_data = get(data_all, id.key, nothing)
                if isnothing(previously_added_data)
                    data_all[id.key] = data
                else
                    dim = is_model_data ? :member : :source
                    prev_models = collect(dims(previously_added_data, dim))
                    new_models = collect(dims(data, dim))
                    joint_data = cat(
                        previously_added_data, data;
                        dims=Dim{dim}(vcat(Array(prev_models), Array(new_models)))
                    )
                    joint_meta = joinMetadata(
                        previously_added_data.metadata, 
                        data.metadata,
                        is_model_data
                    )
                    data_all[id.key] = rebuild(joint_data; metadata = joint_meta)
                    push!(ids_all, id)
                end
            end
        end
    end
    result =  Data(base_path = base_path, ids = ids_all, data = data_all)
    @info "The following data was found and loaded: " result.ids
    if is_model_data && common_models_across_vars
        @info "only retain models shared across all variables"
        result = getCommonModelsAcrossVars(result, :member)
    end
    if length(result.data) != length(result.ids)
        @warn "length of data and ids doesnt match!"
    end
    return result
end
    

"""
    getCommonModelsAcrossVars(modelData::Data, dim::Symbol)

Return only those models (on level of ensemble members) for which there is 
data for all variables.

# Arguments
- `modelData`:
"""
function getCommonModelsAcrossVars(modelData::Data, dim::Symbol)
    data_all = deepcopy(modelData.data)
    ids_all = copy(modelData.ids)
    variables = unique(map(id -> id.variable, ids_all))
    shared_models =  nothing
    for var in variables
        modelDict = filter(((k,v),)-> occursin(var, k), data_all)
        # iterate over all combinations (of diagnostics/statistics) with current variable
        for (_, data_var) in modelDict
            models = collect(dims(data_var, dim))
            if isnothing(shared_models)
                shared_models = models
            else
                shared_models = intersect(shared_models, models)
            end
        end
    end
    for id in map(id -> id.key, ids_all)
        # TODO: only subset if it is necessary at all 
        data_all[id] = subsetModelData(data_all[id], Array(shared_models));
    end
    result = Data(
        base_path = modelData.base_path, ids = modelData.ids, data = data_all
    )
    #alignIDsFilteredData!(result)
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
    loadWeightsAsDimArray(path_to_file::String, key_weights::String)

# Arguments:
- `data`: NCDataset containing weights, which have a single dimension
- `key_weights`: name of weights to load; 'wP' (performance weights), 'wI'
(independence weights), 'w' (overall weights)
"""
function loadWeightsAsDimArray(data::NCDataset, key_weights::String)
    src_name = dimnames(data[key_weights])[1]
    sources = Array(data[src_name])
    arr = DimArray(
        Array(data[key_weights]), 
        (Dim{Symbol(src_name)}(sources)), metadata = Dict(data.attrib)
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