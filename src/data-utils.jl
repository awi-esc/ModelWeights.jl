# these are identical across models
GLOBAL_METADATA_KEYS = [
    "realm",
    "variable_id",
    "experiment_id",
    "product",
    "software",
    "mip_era",
    "parent_mip_era",
    "external_variables",
    "project_id",
    "Conventions",
    "activity_id",
    "sub_experiment_id",

    "standard_name",
    "coordinates",
    "_FillValue",
    "units",
    "long_name",

    "grid",
    "frequency",
    "cell_methods"
    ];
    
# these differ across models and will be saved as list
LOCAL_METADATA_KEYS = [
    "source_id",
    "institution_id",
    "variant_label",
    "grid_label",
];

function checkMetadata(data)
    filtered = filter(((k,v),)->isa(v, Vector) && length(v) != length(dims(data, :model)) , data.metadata);
    if !isempty(filtered) && any(x -> length(x) != 0, values(filtered))
        @warn "Check Metadata, there are fields which were supposed to be identical across models, but werent!" filtered
    end
end

function hasFlawedMetadata(attributes, filename)
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


function createMetaDict(list_keys::Vector{String} = LOCAL_METADATA_KEYS)
    meta = Dict();
    for k in list_keys
        meta[k] = [];
    end
    return meta
end


function verifyMetadata!(meta::Dict, model_dim::Number)
    indices = findall(x -> x != model_dim, map(length, values(meta)));
    for i in indices 
        @warn "The following metadata key seems to be missing for some files: " collect(keys(meta))[i];
    end
    vals = map(k -> (k, unique(meta[k])), collect(keys(meta)))
    unequal_values = filter(
        ((k, l), ) -> length(l) > 1 && !(k in SimilarityWeights.LOCAL_METADATA_KEYS), vals
    )
    if length(unequal_values) > 0
        msg = join(map(x -> x[1], unequal_values), ", ", ", ") 
        @warn "values across models that aren't identical: " msg
    end
    for (k, l) in vals 
        if !(k in SimilarityWeights.LOCAL_METADATA_KEYS)
            meta[k] = l[1]
        end
    end
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
        meta = createMetaDict();
        for (root, dirs, files) in walkdir(pathToData; follow_symlinks=true)
            ncFiles = filter(x->endswith(x, ".nc"), files);
            for file in ncFiles
                addFile = true;
                # only include files that contain all names given in 'included'
                if length(included) != 0
                    if !all([occursin(name, file) for name in included])
                        addFile = false;
                    end
                end
                if addFile
                    filename = split(basename(file), ".nc")[1];
                    ds = NCDataset(joinpath(root, file));
                    dsVar = ds[climVar];

                    attributes = deepcopy(ds.attrib);
                    attributes = merge(Dict(dsVar.attrib), attributes);
                    if hasFlawedMetadata(attributes, filename)
                        continue
                    end
                    attributes = filter(((k,v),) -> k in GLOBAL_METADATA_KEYS || k in LOCAL_METADATA_KEYS, attributes);
                    meta = mergewith(appendValuesDicts, meta, Dict(attributes));
                    if "model_id" in keys(ds.attrib)
                        push!(sources, ds.attrib["model_id"]);
                    else
                        push!(sources, filename)#split(filename, "_")[2])
                    end
                    if climVar == "amoc"
                        if "season_number" in keys(ds.dim)
                            dim1 = Dim{:season}(collect(dsVar["season_number"][:]));
                            push!(data, DimArray(Array(dsVar), (dim1)));
                        else
                            push!(data, DimArray(Array(dsVar), ()));
                        end
                    else
                        # TODO: hur has three dimensions, for now just lon-lat dimensions supported
                        if length(size(dsVar)) != 2
                            throw(ArgumentError(join(["only variables that have ONLY dimensions lon, lat are supported,", 
                                                dsVar.attrib["standard_name"], "has size", size(dsVar)], " ", " ")))
                        end
                        dim1 = Dim{:lon}(collect(dsVar["lon"][:]));
                        dim2 = Dim{:lat}(collect(dsVar["lat"][:]));          
                        push!(data, DimArray(Array(dsVar), (dim1, dim2)));
                    end
                end
            end
        end
        dimData = cat(data...; dims = (Dim{:model}(sources)));
        
        verifyMetadata!(meta, length(sources));
        dimData = rebuild(dimData; metadata = meta);

        checkMetadata(dimData);
        dataAllVars[climVar] = dimData;
    end


    return dataAllVars
end

function appendValuesDicts(val1, val2)
    if isa(val1, Vector) && isa(val2, Vector) 
        return vcat(val1, val2)
    elseif isa(val1, Vector)
        return push!(val1, val2)
    elseif isa(val2, Vector)
        return push!(val2, val1)
    else
        return [val1, val2]
        # if val1 == val2
        #     return val1
        # else
        #     return [val1, val2]
        # end
    end
end