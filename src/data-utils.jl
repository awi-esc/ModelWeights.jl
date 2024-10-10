import YAML
using DimensionalData
using Interpolations

@kwdef mutable struct Config 
    base_path::String
    path_to_recipes::String
    prefix_recipes::String
    target_dir::String
    diagnostics::Vector{String}
    variables::Vector{String}
    name_ref_period::String
    name_full_period::String
    name_obs_period::String
    models_project_name::String
    obs_data_name::String
    weights_variables::Union{Nothing, Dict{String, Dict{String,Number}}}
    weight_contributions::Union{Nothing, Dict{String, Number}}
end


const MODEL_MEMBER_DELIM = "#"

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
        # add mapping from model (ensemble) names to indices in metadata arrays
        meta["ensemble_names"] = Array(dims(data, :model))
        meta["indices_map"] = getIndicesMapping(meta["ensemble_names"]);
    end

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
    loadPreprocData(diagnosticVarToPaths::Dict{String, Dict{String, String}}, dataIncluded=[])

Loads all preprocessed (e.g., by ESMValTool) .nc files (each model is a different .nc file) for each climate climate variable
inside the respective subdirectory (e.g. tos) into a single DimArray with dimensions lon (longitude), lat (latitude) and model. 

If length(dataIncluded) != 0 only those .nc files are considered whose filenames contain all elements of dataIncluded, 
e.g. if dataIncluded=['ERA5'] only ERA5 data will be included (files with ERA5 in their filename).

# Arguments
-`diagnosticVarToPaths`: dictionary mapping from diagnostics (e.g., CLIM) to a
dictionary from climate variables as short names (e.g., tas) to path where 
respective (preprocessed) data of that climate variable is stored.
-`included`: Array that contains Strings that must occur in filenames of loaded
data. If only a certain model should be loaded, this is specified here, e.g.
['AWI-ESM-1-1-LR'] would load only data from this particular model.
"""
function loadPreprocData(
    diagnosticVarToPaths::Dict{String, Dict{String, String}},
    included::Vector{String}=[]
)
    dataAllVars = Dict{String, Dict{String, DimArray}}();
    diagnostics = keys(diagnosticVarToPaths)
    for diagnostic in diagnostics
        varToPaths = diagnosticVarToPaths[diagnostic];
        dataAllVars[diagnostic] = Dict{String, DimArray}();
        for climVar in keys(varToPaths)
            pathToData = varToPaths[climVar]
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
                    else # observational data does not have model name
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
                    # end
                    # print("i=" * string(i) * ": ")
                    # println(string(length(get(meta, "initialization_index", []))))
                end
            end
            dimData = cat(data...; dims = (Dim{:model}(sources)));
            updateMetadata!(meta, dimData);
            dimData = rebuild(dimData; metadata = meta);

            # set model names in model dimension to full model name which was 
            # built in updateMetadata! (cannot be done before since missing 
            # files etc. have to be accounted for first)

            if !isempty(get(meta, "full_model_names", []))
                dimData = set(dimData, :model => meta["full_model_names"])
            end
            dataAllVars[diagnostic][climVar] = dimData;
        end
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
    buildPathsToVarData(config::Config, name_time_period::String)

Returns a mapping from diagnostic (e.g. CLIM) to variable (e.g. tas) to the
paths where the data is stored. The assumed structure of the data is this: 
config.base_path contains the path to the directory that contains a directory 
called 'preproc', which contains a directory for every time period considered 
(e.g. historical1, historical) which in turn contains a directory for every 
combination of variable and diagnostic (e.g., tas_CLIM).
"""
function buildPathsToVarData(config::Config, name_time_period::String)
    var_to_path = Dict{String, Dict{String, String}}();
    for diagnostic in config.diagnostics
        var_to_path[diagnostic] = Dict();
        for var in config.variables
            var_diagnostic = var * "_" * diagnostic
            var_to_path[diagnostic][var] = joinpath(
                config.base_path,
                "preproc",
                name_time_period,
                var_diagnostic
            )
        end
    end
    return var_to_path
end


function validateConfig(path_config::String)
    config_yaml = YAML.load_file(path_config);
    weights_variables = get(config_yaml, "weights_variables", nothing);
    if !isnothing(weights_variables)
        weights_variables = convert(
            Dict{String, Dict{String, Number}}, 
            weights_variables
        )
    end
    # one of base_path and path_to_recipes must be provided
    base_path = get(config_yaml, "base_path", "")
    path_to_recipes = get(config_yaml, "path_to_recipes", "")
    if isempty(base_path) && isempty(path_to_recipes)
        msg = "Either 'base_path' (path to directory that contains preproc-directory) or 'path_to_recipes' (path to directory that contains for every variable a subfolder which in turn contain the preproc-directories.)";
        throw(ArgumentError(msg))
    end
    config = Config(
        base_path = base_path,
        path_to_recipes = path_to_recipes,
        prefix_recipes = get(config_yaml, "prefix_recipes", ""),
        target_dir = joinpath(config_yaml["target_dir"], getCurrentTime()),
        
        diagnostics = config_yaml["diagnostics"],
        variables = config_yaml["variables"],
        
        name_ref_period = get(config_yaml, "name_ref_period", ""),
        name_full_period = get(config_yaml, "name_full_period", ""),
        name_obs_period = get(config_yaml, "name_obs_period", ""),
        
        models_project_name = config_yaml["models_project_name"],
        obs_data_name = get(config_yaml, "obs_data_name", ""),
        
        weights_variables = weights_variables,
        weight_contributions = get(config_yaml, "weight_contributions", nothing)
    ) 
    # TODO: add checks consistency, paths for all specified variables there and exist, etc.
    return config
end

"""
    getCommonModelsAcrossVars(modelData::Dict{String, Dict{String, DimArray}})

Keep only those models for which there is data for all variables.

# Arguments
- `modelData`: maps from diagnostic to climate variable to respective data
"""
function getCommonModelsAcrossVars(modelData::Dict{String, Dict{String, DimArray}})
    data_all = deepcopy(modelData);
    diagnostics = keys(modelData)
    shared_models =  nothing
    for diagnostic in diagnostics
        for var in keys(modelData[diagnostic])
            data_var = modelData[diagnostic][var];
            meta = data_var.metadata;
            models = getUniqueModelIds(meta, getCMIPModelsKey(meta));
            data_all[diagnostic][var].metadata["full_model_names"] = models;
            if isnothing(shared_models)
                shared_models = models
            else
                shared_models = intersect(shared_models, models)
            end
        end
    end
    for diagnostic in diagnostics
        for var in keys(modelData[diagnostic])
            data_all[diagnostic][var] = keepModelSubset(data_all[diagnostic][var], shared_models);
        end
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
    data.metadata["indices_map"] = getIndicesMapping(data.metadata["ensemble_names"])
    return data
end


"""
    loadDataFromConfig(config::Config, name_time_period::String, name_data::String)

Load climate model data for variables and time periods defined in 'config'. 

# Arguments:
- `config`: Config generated from config file (e.g, output of fn validateConfig)
- `name_time_period`: e.g. 'historical', 'historical1'
- `name_data`: e.g. CMIP6, CMIP5, ERA5
"""
function loadDataFromConfig(config::Config, name_time_period::String, name_data::String) 
    pathsDict = Dict{String, Dict{String, String}}();
    config = deepcopy(config);
    if !isempty(config.path_to_recipes)
        variables = deepcopy(config.variables)
        for var in variables
            dirname = var
            if !isempty(config.prefix_recipes)
                dirname = config.prefix_recipes * "_" * var
            end
            config.base_path = joinpath(config.path_to_recipes, dirname)
            config.variables = [var]
            newDict = buildPathsToVarData(config,  name_time_period)
            for diagnostic in keys(newDict)
                merged_subdic = get(pathsDict, diagnostic, Dict{String, String}())
                pathsDict[diagnostic] = merge(merged_subdic, newDict[diagnostic])
            end
        end
    else
        pathsDict = buildPathsToVarData(config,  name_time_period)
    end
    data = loadPreprocData(pathsDict, [name_data]);
    return data
end

function saveWeights(
        weights::DimVector,
        target_dir::String, 
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
each vector only contains N entries where N is the number of Ensembles/Models (i.e. without the unique ensemble members).
If the metadata for the ensemble members of a Model/Ensemble differ across the members, the respective 
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


function loadWeightsAsDimArray(path_to_file::String)
    data = NCDataset(path_to_file)
    models = Array(data["model"])
    arr = DimArray(Array(data["weight"]), (Dim{:model}(models)), metadata = Dict(data.attrib))
    return arr
end


"""
    loadModelData(config::Config)

Load model data for full and/or reference period as defined in config.
Model data is only included if the model has data for all variables. 
If full and reference period are loaded, only the data of the models that 
they have in common are loaded.
"""
function loadModelData(config::Config)
    modelDataFull = nothing
    modelDataRef = nothing
    load_full_period = !isempty(config.name_full_period);
    load_ref_period = !isempty(config.name_ref_period);
    
    if load_full_period
        modelDataFull = loadDataFromConfig(config, config.name_full_period, config.models_project_name);
        modelDataFull = getCommonModelsAcrossVars(modelDataFull);
    end
    if load_ref_period
        modelDataRef = loadDataFromConfig(config, config.name_ref_period, config.models_project_name);
        modelDataRef = getCommonModelsAcrossVars(modelDataRef);
    end

    if load_full_period && load_ref_period
        modelDataFull, modelDataRef = getSharedModelData(modelDataFull, modelDataRef);
    end

    return (modelDataFull, modelDataRef)
end

"""
    computeInterpolatedWeightedQuantiles

This implementation follows the one used by Brunner et al.
"""
function computeInterpolatedWeightedQuantiles(quantiles, vals, weights=nothing)
    if isnothing(weights)
        weights = ones(length(vals));
    end
    indicesSorted = Array(sortperm(vals)); # gives indices of smallest to largest data point
    weightsSorted = weights[indicesSorted];
    weightedQuantiles = cumsum(weightsSorted) - 0.5 * weightsSorted;
    weightedQuantiles = reshape(weightedQuantiles, length(weightedQuantiles), 1);
    weightedQuantiles = (weightedQuantiles .- minimum(weightedQuantiles)) ./ maximum(weightedQuantiles);
    
    interp_linear = Interpolations.linear_interpolation(
        vec(weightedQuantiles), 
        vals[indicesSorted],
        extrapolation_bc=Interpolations.Line()
    );
     
    return interp_linear(quantiles)
end



"""
    getUncertaintyRanges(data::DimArray, weights::DimArray, quantiles::Vector{Number})
    
# Arguments:
- `data`: has dimensions 'time', 'model'
- `weights`: has dimension 'model', sum up to 1
- `quantiles`: Vector with two entries (btw. 0 and 1) [lower_bound, upper_bound]
"""
function getUncertaintyRanges(
    data::DimArray,
    w::DimArray,
    quantiles=[0.167, 0.833]
)
    unweightedRanges = [];
    weightedRanges = [];
    for t in dims(data, :time)
        lower, upper = computeInterpolatedWeightedQuantiles(
            quantiles, Array(data[time = At(t)]), w
        );
        push!(weightedRanges, [lower, upper]);
        lower, upper = computeInterpolatedWeightedQuantiles(
            quantiles, Array(data[time = At(t)])
        );
        push!(unweightedRanges, [lower, upper]);
    end

    return (weighted=weightedRanges, unweighted=unweightedRanges)
end

"""
    renameModelDimsFromMemberToEnsemble(data::DimArray, dim_names::Vector{String})

# Arguments:
- `data`
- `dim_names`: names of dimensions to be changed, e.g. 'model', 'model1', etc.
"""
function renameModelDimsFromMemberToEnsemble(data::DimArray, dim_names::Vector{String})
    for dim in dim_names
        unique_members = dims(data, Symbol(dim))
        ensembles = map(x -> split(x, MODEL_MEMBER_DELIM)[1], unique_members)
        data = set(data, Symbol(dim) => ensembles)
    end
    return data
end