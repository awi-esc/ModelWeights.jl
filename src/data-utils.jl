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


function debugMetadata(data)
    filtered = filter(((k,v),)->isa(v, Vector) , data.metadata);
    if !isempty(filtered)
        @debug "Values of metadata differ for these fields:" keys(filtered)
    end
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


function getUniqueModelIds(meta::Dict{String, Union{String, Array}}, key_model_names::String)
    mip_era = get(meta, "mip_era", get(meta, "project_id", ""))
    models = [];
    if mip_era == "CMIP5"
        # TODO add
        @warn "For CMIP5, full model names are not generated in metadata."
    else
        names = meta[key_model_names];
        variants = meta["variant_label"];
        grids = meta["grid_label"]; 
        models = map(
            x->join(x, "__", "__"), 
            zip(names, variants, grids)
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
    meta::Dict{String, Union{Array, String}}, 
    data::DimArray, 
    indices_ignored::Array
)
    model_dim = length(dims(data, :model))
    for (i, key) in enumerate(keys(meta))
        values = meta[key]; 
        # add missing values for the last added files 
        n =  model_dim - length(values)
        for i in range(1, n)
            push!(values, missing)
        end
        #println("i: " * string(i) * " " * string(length(meta[key])))
        deleteat!(values, indices_ignored)
        # if none was missing and all have the same value, just use a string
        if !any(ismissing, values) && length(unique(values)) == 1
            meta[key] = string(values[1])
        end
    end
    mip_era = get(meta, "mip_era", get(meta, "project_id", ""))
    if !isempty(mip_era)
        key_model_name = ifelse(mip_era == "CMIP6", "source_id", "model_id")
        meta["full_model_names"] = getUniqueModelIds(meta, key_model_name)
    end
end


function getMetadataCombinedModels(metadata::Dict, climVar::String)
    models_key = getCMIPModelsKey(metadata);
    meta_shared = filter(((k,v),)->!isa(v, Vector), metadata);
    get!(meta_shared, "variables", climVar);
    get!(meta_shared, "full_model_names", metadata["full_model_names"]);
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
        meta = Dict{String, Union{String, Array}}();
        files = readdir(pathToData; join=true);
        ncFiles = filter(x->endswith(x, ".nc"), files);
        indices_ignored = [];
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
                if warnIfFlawedMetadata(attributes, filename)
                    push!(indices_ignored, i)
                    continue
                end

                # update metadata
                for key in keys(attributes)
                    values = get!(meta, key, []);
                    # fill up vector
                    n = i - length(values) - 1;
                    for j in range(1, n)
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
                # println(string(length(get(meta, "model_doi_url", []))))
            end
        end
        dimData = cat(data...; dims = (Dim{:model}(sources))); 
        updateMetadata!(meta, dimData, indices_ignored);
        dimData = rebuild(dimData; metadata = meta);
        # TODO: if log-level set to debug: debugMetadata(dimData);
        dataAllVars[climVar] = dimData;
    end
    return dataAllVars
end

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
    else
        return [val1, val2]
    end
end


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
        models = getUniqueModelIds(meta, getCMIPModelsKey(meta));
        data_all[var].metadata["full_model_names"] = models;
        if isnothing(shared_models)
            shared_models = models
        else
            shared_models = intersect(shared_models, models)
        end
    end
    keepModelSubset!(data_all, shared_models);
    return data_all
end


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

function keepModelSubset!(data::Dict{String, DimArray}, shared_models::Vector{String})
    for var in keys(data);
        data_var = data[var]
        indices = findall(m -> m in shared_models, data_var.metadata["full_model_names"]);
        data_var = data_var[model = indices];
        meta = data_var.metadata;
        # adapt metadata accordingly
        meta["full_model_names"] = meta["full_model_names"][indices];
        models_key = getCMIPModelsKey(meta)
        meta[models_key] = meta[models_key][indices];
        attributes = filter(x -> !(x isa String), keys(meta));
        for key in attributes
            meta[key] = meta[key][indices];
        end
        data[var] = data_var;
    end
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
