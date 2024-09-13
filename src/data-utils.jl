import YAML
using DimensionalData

@kwdef struct Config 
    base_path::String
    target_dir::String
    experiment::String
    prefix_var_folders::String = ""
    variables::Vector{String}
    name_ref_period::String
    name_full_period::String
    models_project_name::String
    obs_data_name::String
    weights_variables::Dict{String, Dict{String,Number}}
    weight_contributions::Dict{String, Number}
end


function warnIfFlawedMetadata(attributes, filename)
    isFlawed = false;
    if "branch_time_in_parent" in keys(attributes) && isa(attributes["branch_time_in_parent"], String)
        @warn "Branch_time_in_parent is a string, excluded file:" filename
        isFlawed = true
    elseif "branch_time_in_child" in keys(attributes) && isa(attributes["branch_time_in_child"], String)
        @warn "Branch_time_in_child is a string, excluded file:" filename
        isFlawed = true
    end
    return isFlawed
end


function buildCMIP5EnsembleMember(realizations, initializations, physics, n)
    function concat(elem, prefix)
        if !(elem isa Vector)
            return [prefix * string(elem) for _ in range(1, n)]
        else 
            return map(x -> prefix * string(x), elem)
        end
    end
    
    rips = [concat(realizations, "r"), concat(initializations, "i"), concat(physics, "p")];
    variants = [join(k, "") for k in zip(rips...)]
    return variants
end

function getUniqueModelIds(
    meta::Dict{String, Union{String, Array, Dict}},
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
            x->join(x, "__"), 
            zip(model_names, variants)
        );
        @warn "For CMIP5, full model names dont include grid."
    else
        variants = meta["variant_label"];
        grids = meta["grid_label"]; 
        models = map(
            x->join(x, "__", "__"), 
            zip(model_names, variants, grids)
        );
    end
    return models
end


"""
    updateMetadata!(
        meta::Dict{String, Union{Array, String}}, 
        data::DimArray, 
        indices_ignored::Array
    )

Update metadata 'meta' s.t. data of ignored files is removed and attributes
that were only present in some files/models are set to missing. Further the key 
'full_model_names' is added which contains for every file/model the unique 
identifier consisting of variant_label, grid_label and model_name. 
"""
function updateMetadata!(
    meta::Dict{String, Union{Array, String, Dict}}, 
    data::DimArray
)
    model_dim = length(dims(data, :model))
    for key in keys(meta)
        values = meta[key]; 
        # add missing values for the last added files 
        n =  model_dim - length(values)
        for _ in range(1, n)
            push!(values, missing)
        end
        #println("i: " * string(i) * " " * string(length(meta[key])))
        # if none was missing and all have the same value, just use a string
        if !any(ismissing, values) && length(unique(values)) == 1
            meta[key] = string(values[1])
        end
    end
    # for model data only
    if !isempty(get(meta, "source_id", get(meta, "model_id", "")))
        key_model_name = getCMIPModelsKey(meta);
        meta["full_model_names"] = getUniqueModelIds(meta, key_model_name)
    end
    # add mapping from model names to indices in metadata arrays
    meta["indices_map"] = getIndicesMapping(Array(dims(data, :model)));
end


function getIndicesMapping(names)
    mapping = Dict();
    for name in names
        mapping[name] = findall(x -> x == name, names)
    end
    return mapping
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
    #meta_shared = filter(((k,v),)->!isa(v, Vector), metadata);
    models_key = getCMIPModelsKey(metadata);
    meta_shared = filter(((k,v),) -> isa(v, String), metadata);
    meta_shared["full_model_names"] = deepcopy(metadata)["full_model_names"];
    meta_shared[models_key] = unique(metadata[models_key]);
    return meta_shared
end


"""
    loadPreprocData(climateVarsToPaths::Dict{String,String}, dataIncluded=[])

Loads all preprocessed (e.g., by ESMValTool) .nc files (each model is a different .nc file) for each climate climate variable
inside the respective subdirectory (e.g. tos) into a single DimArray with dimensions lon (longitude), lat (latitude) and model. 

If length(dataIncluded) != 0 only those .nc files are considered whose filenames contain all elements of dataIncluded, 
e.g. if dataIncluded=['ERA5'] only ERA5 data will be included (files with ERA5 in their filename).

# Arguments
'climateVarsToPaths': dictionary mapping from climate variables as short names (e.g., tas) to path where respective (preprocessed) data is stored
'included': Array that contains Strings that must occur in filenames of loaded data. If only a certain model should be loaded, this is specified here, e.g. 
            ['AWI-ESM-1-1-LR'] would load only data from this particular model.
"""
function loadPreprocData(climVarsToPaths::Dict{String, String}, included::Vector{String}=[])
    
    dataAllVars = Dict{String, DimArray}();
    for climVar in keys(climVarsToPaths)
        pathToData = climVarsToPaths[climVar];
        if !isdir(pathToData)
            throw(ArgumentError(pathToData * " does not exist!"))
        end
        data = [];
        sources = [];
        meta = Dict{String, Union{String, Array, Dict}}();
        files = readdir(pathToData; join=true);
        ncFiles = filter(x->endswith(x, ".nc"), files);
        nbIgnored = 0;
        for (i, file) in enumerate(ncFiles)
            addFile = true;
            # only include files that contain all names given in 'included'
            if length(included) != 0
                if !all([occursin(name, file) for name in included])
                    addFile = false;
                end
            end
            if addFile
                filename = split(basename(file), ".nc")[end-1];
                ds = NCDataset(file);
                dsVar = ds[climVar];

                attributes = merge(Dict(deepcopy(dsVar.attrib)), Dict(deepcopy(ds.attrib)));
                if climVar == "msftmz"
                    attributes = merge(attributes, Dict(deepcopy(ds["sector"].attrib)));
                end
                if SimilarityWeights.warnIfFlawedMetadata(attributes, filename)
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
                else
                    # for observational data
                    push!(sources, filename)
                end

                # if climVar == "amoc"
                #     if "season_number" in keys(ds.dim)
                #         dim1 = Dim{:season}(collect(dsVar["season_number"][:]));
                #         push!(data, DimArray(Array(dsVar), (dim1)));
                #     else
                #         push!(data, DimArray(Array(dsVar), ()));
                #     end
                # else
                dimension_names = dimnames(dsVar)
                dimensions = []
                for d in dimension_names
                    if d in ["bnds", "string21"]
                        continue
                    end
                    push!(dimensions, Dim{Symbol(d)}(collect(dsVar[d][:])))
                end      
                push!(data, DimArray(Array(dsVar), Tuple(dimensions)));
                # end
                # print("i=" * string(i) * ": ")
                # println(string(length(get(meta, "initialization_index", []))))
            end
        end
        dimData = cat(data...; dims = (Dim{:model}(sources))); 
        updateMetadata!(meta, dimData);
        dimData = rebuild(dimData; metadata = meta);
        dataAllVars[climVar] = dimData;
    end
    return dataAllVars
end

"""
    appendValuesDicts(val1, val2)

Combine values of two different dictionaries iteratively. If both values are 
Vectors, they should be identical and only one of them is added. If they
aren't identical, a warning is triggered and both vectors are concatenated 
into one big Vector.

# Arguments:
- `val1`: a Vector, a Dictionary or a single value (e.g. String, Number, etc.)
- `val2`: a Vector, a Dictionary or a single value (e.g. String, Number, etc.)
"""
function appendValuesDicts(val1, val2)
    if isa(val1, Vector) && isa(val2, Vector) & !(isempty(val1) || isempty(val2))
        if val1 != val2 
            @warn "Two arrays merged that weren't identical! (usuallly in metadata)"
            @warn val1
            @warn val2
            return vcat(val1, val2)
        else
            return val1 
        end
    elseif isa(val1, Vector)
        return push!(val1, val2)
    elseif isa(val2, Vector)
        return push!(val2, val1)

    elseif isa(val1, Dict) && isa(val2, Dict)
        if val1 != val2
            @warn "Two different dictionaries (check indices_map)"
            @warn val1 
            @warn val2
            return vcat(val1, val2)
        else
            return val1
        end

    elseif isa(val1, Dict) || isa(val2, Dict)
        @warn val1 
        @warn val2
        throw(ArgumentError("Dictionary merged with something else!"))
    else
        return [val1, val2]
    end
end

"""
    buildPathsToVarData(config::Config, period::String)
"""
function buildPathsToVarData(config::Config, period::String)
    prefix = config.prefix_var_folders;
    var_to_path = Dict{String, String}();
    for var in config.variables
        folder = ifelse(isempty(prefix), var, prefix * "_" * var);
        path_data = joinpath(
            config.base_path,
            config.experiment, 
            folder,
            "preproc"
        )
        var_to_path[var] = joinpath(
            path_data, 
            "climatology_" * period,
            var
        )
    end
    return var_to_path
end


function validateConfig(path_config::String)
    config_yaml = YAML.load_file(path_config);
    config = Config(
        base_path = config_yaml["base_path"],
        target_dir = joinpath(config_yaml["target_dir"], getCurrentTime()),
        experiment = config_yaml["experiment"],
        prefix_var_folders = config_yaml["prefix_var_folders"],
        variables = config_yaml["variables"],
        name_ref_period = config_yaml["name_ref_period"],
        name_full_period = config_yaml["name_full_period"],
        models_project_name = config_yaml["models_project_name"],
        obs_data_name = config_yaml["obs_data_name"],
        weights_variables = convert(
            Dict{String, Dict{String, Number}}, 
            config_yaml["weights_variables"]
        ), 
        weight_contributions = config_yaml["weight_contributions"]
    )
    # TODO: add checks consistency, paths for all specified variables there and exist, etc.
    return config
end

"""
    getCommonModelsAcrossVars(modelData::Dict{String, DimArray})

Keep only those models for which there is data for all variables.

# Arguments
- `modelData`: maps from climate variable to respective data
"""
function getCommonModelsAcrossVars(modelData::Dict{String, DimArray})
    data_all = deepcopy(modelData);
    variables = keys(modelData);
    
    shared_models =  nothing
    for var in variables
        data_var = modelData[var];
        meta = data_var.metadata;
        models = SimilarityWeights.getUniqueModelIds(meta, SimilarityWeights.getCMIPModelsKey(meta));
        data_all[var].metadata["full_model_names"] = models;
        if isnothing(shared_models)
            shared_models = models
        else
            shared_models = intersect(shared_models, models)
        end
    end
    for var in variables
        data_all[var] = SimilarityWeights.keepModelSubset(data_all[var], shared_models);
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
    keepMetadataSubset!(meta::Dict, indices::Vector{Number})

# Arguments:
- `meta::Dict`: metadata dictionary
- `indices::Vector{Number}`: indices of data to be kept
"""
function keepMetadataSubset!(meta::Dict, indices::Vector{Int64})
    attributes = filter(x -> meta[x] isa Vector, keys(meta));
    for key in attributes
        meta[key] = meta[key][indices];
    end
end

"""
    keepModelSubset(data::Dict{String, DimArray}, shared_models::Vector{String})

Retain data only from models in `shared_models`. Takes care of metadata.
"""
function keepModelSubset(data::DimArray, shared_models::Vector{String})
    indices = findall(m -> m in shared_models, data.metadata["full_model_names"]);
    @assert length(indices) == length(shared_models)
    keepMetadataSubset!(data.metadata, indices);
    data = data[model = indices];
    data.metadata["indices_map"] = getIndicesMapping(Array(dims(data, :model)));
    return data
end


"""
    loadModelDataFromConfig(config::Config, key_period::String)

Load climate model data for variables and time periods defined in 'config'. 

# Arguments:
- `config`: Config generated from config file (e.g, output of fn validateConfig)
- `key_period`: e.g. 'name_ref_period', 'name_full_period'
"""
function loadDataFromConfig(config::Config, key_period::String, key_data_name::String)
    pathsDict = buildPathsToVarData(config,  getproperty(config, Symbol(key_period)));
    data = loadPreprocData(pathsDict, [getproperty(config, Symbol(key_data_name))]);
    return data
end


function saveWeights(
        weights::DimVector,
        averages::Dict{String, Dict{String, DimArray}},
        target_dir::String, 
        target_fn::String="output.nc"
    )
    if !isdir(target_dir)
        mkpath(target_dir)
    end
    path_to_target = joinpath(target_dir, target_fn);
    # TODO: only works if not yet created!
    ds = NCDataset(path_to_target, "c")

    defDim(ds, "model", length(dims(weights, :model)))
    # global attributes
    for (k, v) in weights.metadata
        ds.attrib[k] = v
    end
    v = defVar(ds, "weight", Float64, ("model",))
    v[:] = Array(weights)
    close(ds)

    # TODO: adapt for non lon,lat -data!
    # TODO: also add/log computed averages (weighted/uweighted)
    # for avg in ["weighted", "unweighted"]
    #     for (climVar, values) in averages[avg]
    #         name = join([avg, "avg", climVar], "_", "_")
    #         var = defVar(ds, name, Float32, ("lon", "lat"))
    #         # Replace missing with the fill value
    #         attput(var, "_FillValue", 1.0e20) 
    #         var[:,:] = coalesce.(Array(values), missing)  
    #     end
    # end

    @info "saved data to " path_to_target
end

"""
    updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)

Vectors in metadata 'meta' have to refer to different models. These are now summarized such that
each vector only contains N entries where N is the number of Models (without ensemble members).
If the metadata for the ensemble members of a Model differ across the members, the respective 
entry in the vector will be a vector itself. 
"""
function updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)
    meta_new = filter(((k,v),) -> !(v isa Vector), meta);    
    attributes = filter(x -> meta[x] isa Vector, keys(meta));
    attribs_diff_across_members = [];
    for key in attributes
        for model in dims(grouped_data, :model)
            indices = meta["indices_map"][model]
            vals = get!(meta_new, key, [])       
            val_ensemble = unique(meta[key][indices]);
            if length(val_ensemble) != 1
                push!(vals, val_ensemble)
                push!(attribs_diff_across_members, key)
            else
                push!(vals, val_ensemble[1])
            end
        end
    end
    if !isempty(attribs_diff_across_members)
        @warn "metadata attributes that differ across ensemble members (ok for some!)" unique(attribs_diff_across_members)
        @assert !(getCMIPModelsKey(meta) in attribs_diff_across_members)
    end
    return meta_new
end
