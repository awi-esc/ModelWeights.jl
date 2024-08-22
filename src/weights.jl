using NCDatasets
using DimensionalData
using Statistics
using LinearAlgebra

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
    "frequency",

    "standard_name",
    "coordinates",
    "_FillValue",
    "units",
    "long_name"
    ];

# these differ across models and will be saved as list
LOCAL_METADATA_KEYS = [
    "source_id",
    "institution_id",
    "variant_label",
    "cell_methods"
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
        data = [];
        sources = [];
        meta = SimilarityWeights.createMetaDict();
        for (root, dirs, files) in walkdir(pathToData; follow_symlinks=true)
            ncFiles = filter(x->endswith(x, ".nc"), files);
            for file in ncFiles
                addFile = true;
                # if not empty, only includes files that contain all names given in 'included'
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
                    if SimilarityWeights.hasFlawedMetadata(attributes, filename)
                        continue
                    end
                    attributes = filter(((k,v),) -> k in GLOBAL_METADATA_KEYS || k in LOCAL_METADATA_KEYS, attributes);
                    meta = mergewith(SimilarityWeights.appendValuesDicts, meta, Dict(attributes));
                    if "model_id" in keys(ds.attrib)
                        push!(sources, ds.attrib["model_id"]);
                    else
                        push!(sources, split(filename, "_")[2])
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
        dimData = rebuild(dimData; metadata = meta);

        SimilarityWeights.checkMetadata(dimData);
        dataAllVars[climVar] = dimData;
    end


    return dataAllVars
end

function areaWeightedMSE(m1::DimArray, m2::DimArray, mask::DimArray)
    latitudes = dims(m1, :lat);
    areaWeights = cos.(deg2rad.(latitudes));
    areaWeights = DimArray(areaWeights, Dim{:lat}(Array(latitudes)));

    sqDiff = (m1 .- m2).^2;
    weightedValues = areaWeights' .* sqDiff;

    weightMatrix = repeat(areaWeights', length(DimensionalData.dims(m1, :lon)), 1);  
    weightMatrix = ifelse.(mask .== 1, 0, weightMatrix); 
    normalization = sum(weightMatrix);

    mse = sqrt(sum(skipmissing(weightedValues)./normalization));
    return mse
end

"""
    getModelDistMatrix(modelData::DimArray)

Computes the area weighted root mean squared error between model predictions for each pair of models. 

'modelData' has dimensions: lon, lat, model and contains the data for a single variable.

returns a DimArray of size nxn (which is a symmetrical matrix) where n is the number of models in 'modelData'. 
"""
function getModelDistances(modelData::DimArray)
    # Make sure to use a copy of the data, otherwise, it will be modified by applying the mask!!
    data = deepcopy(modelData);
    # only take values where none (!) of the models has infinite values!! (Not just the two that are compared to one another)
    nbModels = length(dims(data, :model));
    maskMissing = dropdims(any(ismissing, data, dims=:model), dims=:model);

    matrixS = zeros(nbModels, nbModels);    
    for (i, model_i) in enumerate(eachslice(data; dims=:model))
        model_i[maskMissing .== 1] .= 0;
        
        for (j, model_j) in enumerate(eachslice(data[:, :, i+1:end]; dims=:model))
            idx = j + i;
            model_j[maskMissing .== 1] .= 0;
            s_ij = areaWeightedMSE(model_i, model_j, maskMissing);
            matrixS[i, idx] = s_ij
        end
    end
    # make symmetrical matrix
    symDistMatrix = matrixS .+ matrixS';
    dim = Array(dims(modelData, :model));
    return DimArray(symDistMatrix, (Dim{:model1}(dim), Dim{:model2}(dim)))
end


"""
    weightDistMatrix(distMatrix::DimArray{T}, weight::T) where T<:Number

Normalizes the matrix by its median and multiplies each entry by 'weight'.
"""
function normalizeAndWeightDistMatrix(distMatrix::Union{DimVector, DimArray}, weight::T=1) where T<:Real
    distMatrix = distMatrix ./ median(distMatrix); 
    weightedDistances = weight .* distMatrix;
    return weightedDistances
end


"""
    normalizeWeightsVariables(weightsVars::Dict{String, Number})

Modifies input vector 'weights' by normalizing it s.t. its elements sum up to 1.
"""
function normalizeWeightsVariables!(weights::Dict{String, Number})
    total = sum(values(weights))
    for key in keys(weights)
        weights[key] /= total 
    end
end

"""
    getIndependenceWeights(data::Dict{String, DimArray}, weightsVars::Dict{String, Number}=Dict{String, Number}())

Computes an independence weight for each model in 'data'. The weights result from taking the average weighted root mean squared errors between pairs of models for each variable. 

'data' is a Dictionary mapping from the climate variable (e.g. tas) to a DimArray with the corresponding data (lon, lat, model). 
'weightsVars' is a Dictionary mapping from the climate variable to a number which is the weight of how much the respective variable contributes to the computed independenceWeight.

returns a DimArray (model1, model2, variable) with the computed independence weights, seperately for each considered variable.
"""
function getIndependenceWeights(data::Dict{String, DimArray}, weightsVars::Dict{String, Number}=Dict{String, Number}())
    variables = keys(data);
    weights = copy(weightsVars);
    if !isempty(weights)
        normalizeWeightsVariables!(weights);
    else
        weights = Dict(zip(variables, ones(length(variables))));
    end

    weightedDistMatrices = [];
    meta = createMetaDict(GLOBAL_METADATA_KEYS);
    meta_not_shared = Dict();
    for climVar in variables
        metadata = data[climVar].metadata;
        meta_shared = filter(((k,v),)->!isa(v, Vector), metadata);
        meta_not_shared[climVar] = filter(((k,v),)->isa(v, Vector), metadata);

        distances = getModelDistances(data[climVar]);
        weight = ifelse(isempty(weights), 1, weights[climVar]);
        weightedDistances = normalizeAndWeightDistMatrix(distances, weight);
        push!(weightedDistMatrices, weightedDistances);
        meta = mergewith(appendValuesDicts, meta, meta_shared);
    end
    meta = merge(meta, meta_not_shared);
    weightsByVars = cat(weightedDistMatrices..., dims = Dim{:variable}(collect(variables)));
    weightsByVars = rebuild(weightsByVars; metadata = meta);
    return weightsByVars
end

"""
    getModelDataDist(models::DimArray, observations::DimArray)

Computes the distance (area-weighted rms error) between model predictions and observations. 

"""
function getModelDataDist(models::DimArray, observations::DimArray)      
    # Make sure to use a copy of the data, otherwise, it will be modified by applying the mask!!
    models = deepcopy(models);
    observations = deepcopy(observations);
    distances = [];
    model_names = [];
    for (i, model_i) in enumerate(eachslice(models; dims=:model))
        model_name = dims(models, :model)[i]  # Access the name of the current model

        maskNbMissing = (ismissing.(observations) + ismissing.(model_i)) .> 0; # observations or model is missing (or both)
        maskedObs = deepcopy(observations);
        maskedObs = dropdims(ifelse.(maskNbMissing .> 0, 0, maskedObs), dims=:model);

        maskedModel = ifelse.(maskNbMissing .> 0, 0, model_i);
        mse = areaWeightedMSE(maskedModel, maskedObs, maskNbMissing);

        push!(model_names, model_name);
        push!(distances, mse);
    end
    return DimArray(distances, (Dim{:model}(model_names)), metadata=models.metadata)
end


"""
    getPerformanceWeights(modelData::Dict{String, DimArray}, obsData::Dict{String, DimArray}, weightsVars::Dict{String, Number}=Dict{String, Number}())

Computes a performance weight for each model in 'modelData'. The weights result from taking the average weighted root mean squared errors between each model and the observational data.
'mdoelData' and 'obsData' are Dictionaries mapping from the climate variable (e.g. tas) to a DimArray with the corresponding data (lon, lat, model). For the observational data, 
the third dimension (model) just contains a single entry (e.g. ERA5).
'weightsVars' is a Dictionary mapping from the climate variable to a number which is the weight of how much the respective variable contributes to the computed performanceWeight.

returns a DimArray (model1) with the computed performance weights, seperately for each considered variable.
"""
function getPerformanceWeights(modelData::Dict{String, DimArray}, 
                               obsData::Dict{String, DimArray},
                               weightsVars::Dict{String, Number}=Dict{String, Number}()
                               )
    variables = keys(modelData);
    weights = copy(weightsVars);
    if !isempty(weights)
        normalizeWeightsVariables!(weights);
    else
        weights = Dict(zip(variables, ones(length(variables))));
    end
    weightedDistMatrices = [];

    meta = createMetaDict(GLOBAL_METADATA_KEYS);
    meta_not_shared = Dict();
    for climVar in variables
        metadata = modelData[climVar].metadata;
        meta_shared = filter(((k,v),)->!isa(v, Vector), metadata);
        meta_not_shared[climVar] = filter(((k,v),)->isa(v, Vector), metadata);

        distances = getModelDataDist(modelData[climVar], obsData[climVar]);
        weight = ifelse(isempty(weights), 1, weights[climVar]);

        weightedDistances = normalizeAndWeightDistMatrix(distances, weight);
        push!(weightedDistMatrices, weightedDistances);
        meta = mergewith(appendValuesDicts, meta, meta_shared);
    end
    meta = merge(meta, meta_not_shared);
    weightsByVar = cat(weightedDistMatrices..., dims = Dim{:variable}(collect(variables)));
    weightsByVar = rebuild(weightsByVar; metadata = meta);
    return weightsByVar
end

function summarizeWeightsAcrossVars(weightsByVar::DimArray)
    summarizedWeights = reduce(+, weightsByVar, dims=:variable);
    return dropdims(summarizedWeights, dims=:variable)
end


""" 
    averageEnsembleVector!(data::DimArray)

For each model, compute the mean across all ensemble members of that model.

# Arguments:
'data': a DimArray with dimensions 'model' and 'variable'

returns a DimArray with dimensions 'variable', 'model'
"""
function averageEnsembleVector(data::DimArray)
    if !hasdim(data, :variable)
        grouped = groupby(data, :model=>identity);
        averages = map(d -> mean(d, dims=:model), grouped);
        models = Array(dims(averages, :model));
        combined = cat(averages..., dims=(Dim{:model}(models)));
    else
        variables = Array(dims(data, :variable));
        results = [];
        for climVar in variables
            # grouped = groupby(data[variable=At(climVar)], :model=>identity);
            grouped = groupby(data[variable = Where(x -> x == climVar)], :model=>identity);
            averages = map(d->mean(d, dims=:model), grouped);
            models = Array(dims(averages, :model));
            avg = cat(averages..., dims=(Dim{:model}(models)));
            push!(results, avg);
        end
        models = unique(Array(dims(data, :model)));
        combined = cat(results...; dims=(Dim{:variable}(variables)));
    end

    return combined
end

"""
    averageEnsembleMatrix(data::DimArray)

Computes the average weights across all members of each model ensemble.

# Arguments:
'data' has dimensions 'model1', 'model2', 'variable'
"""
function averageEnsembleMatrix(data::DimArray)
    if !hasdim(data, :variable)
        grouped = groupby(data, :model1 => identity, :model2 => identity);
        averages = map(d -> mean(mean(d, dims=:model2), dims=:model1)[1,1], grouped);
        averages[LinearAlgebra.diagind(averages)] .= 0;
        combined = averages;
    else
        variables = Array(dims(data, :variable));
        results = [];
        for climVar in variables
            distances = data[variable=Where(x -> x == climVar)];
            grouped = groupby(distances[variable=At(climVar)], :model1 => identity, :model2=>identity);
            averages = map(d-> mean(mean(d, dims=:model2), dims=:model1)[1,1], grouped);
            # Note: set all comparisons of same model to itself to 0, may differ from zero for those models with several ensemble members
            averages[LinearAlgebra.diagind(averages)] .= 0;
            push!(results, averages)
        end
        combined = cat(results..., dims=Dim{:variable}(variables));
    end
    result = rebuild(combined; metadata = data.metadata);
    return result
end

"""
    combineWeights(performanceWeights::DimArray, independenceWeights::DimArray, sigmaD::Number=0.5, sigmaS::Number=0.5)

Combines the RMSEs between pairs of models and the RMSEs between each model and the data into a set of normalized weights, 
one for each model. 

Note that, the given weights are summarized wrt ensemble models, s.t. for each ensemble there is only one value.

# Arguments:
'performanceWeights': DimArray with dimensions 'model', 'variable'
'independenceWeights': DimArray with dimensions 'model1', 'model2' and 'variable'
'sigmaD': Nb. btw. 0 and 1; contribution of performance metric 
'sigmaS': Nb. btw. 0 and 1; contribution of independence metric
"""
function combineWeights(performanceWeights::DimArray, independenceWeights::DimArray, sigmaD::Number, sigmaS::Number)
    wP = averageEnsembleVector(performanceWeights);
    wI = averageEnsembleMatrix(independenceWeights);
    wP = summarizeWeightsAcrossVars(wP);
    wI = summarizeWeightsAcrossVars(wI);
    performance = exp.(-(wP ./ sigmaD).^2);
    # note: (+1 in eq. in paper is from when model is compared to itself since exp(0)=1)
    independence = dropdims(sum(exp.(-(wI ./ sigmaS).^2), dims=:model2), dims=:model2)
    independence = set(independence, :model1 => :model);

    weights = performance ./ independence;
    weightsNormalized = weights ./ sum(weights);
    return weightsNormalized
end

"""
    getWeights(modelData::Dict{String, DimArray}, obsData::Dict{String, DimArray}[, sigmaD::Number=0.5, sigmaS::Number=0.5, 
               weightsPerform::Dict{String, Number}=Dict{String, Number}(), 
               weightsIndep::Dict{String, Number}=Dict{String, Number}()])

    Computes weight for each model in multi-model ensemble according to approach from 
    Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner, Anna L. Merrifield, Ruth Lorenz, and Reto Knutti.
    “Reduced Global Warming from CMIP6 Projections When Weighting Models by Performance and Independence.” 
    Earth System Dynamics 11, no. 4 (November 13, 2020): 995–1012. https://doi.org/10.5194/esd-11-995-2020.

"""
function getWeights(modelData::Dict{String, DimArray}, 
                    obsData::Dict{String, DimArray}, 
                    sigmaD::Number=0.5, 
                    sigmaS::Number=0.5,
                    weightsPerform::Dict{String, Number}=Dict{String, Number}(), 
                    weightsIndep::Dict{String, Number}=Dict{String, Number}()
    )
    wP = getPerformanceWeights(modelData, obsData, weightsPerform);
    wI = getIndependenceWeights(modelData, weightsIndep);
    weights = combineWeights(wP, wI, sigmaD, sigmaS);
    return weights
end

function appendValuesDicts(val1, val2)
    if isa(val1, Vector) && isa(val2, Vector) 
        #print("both are vectors")
        return vcat(val1, val2)
    elseif isa(val1, Vector)
        #print("just arg1 is vector")
        return push!(val1, val2)
    elseif isa(val2, Vector)
        #print("just arg2 is vector")
        return push!(val2, val1)
    else
        #print("none is vector")
        if val1 == val2
            return val1
        else
            return [val1, val2]
        end
    end
end
