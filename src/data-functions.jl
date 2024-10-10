using DimensionalData

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