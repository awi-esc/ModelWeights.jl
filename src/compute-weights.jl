using NCDatasets
using DimensionalData
using Glob
using Statistics

"""
    loadPreprocData(pathToPreprocDir, climateVar, diagnostic=CLIM, dataIncluded=[])

Loads all .nc files (each model is a different .nc file) for variable climateVar and 'diagnostic' inside the respective subdirectory (climateVar_diagnostic)
of the directory located at 'path2PreprocDir' (assuming data had been preprocessed by ESMValTool) into a single DimArray with dimensions lon (longitude), lat (latitude) and model. 

If length(dataIncluded) != 0 only those .nc files are considered whose filenames contain all elements of dataIncluded, 
e.g. if dataIncluded=['ERA5'] only ERA5 data will be included (files with ERA5 in their filename).

'pathToPreprocDir': path to directory with preprocessed data from ESMValTool
'climateVariables': Array of short names of considered variables (e.g. tas, psl)
'diagnostic': the name of the diagnotsic used by ESMValTool
'included': Array that contains Strings that must occur in filenames of loaded data. If only a certain model should be loaded, this is specified here, e.g. 
            ['AWI-ESM-1-1-LR'] would load only data from this particular model.

returns a dictionary from climateVariable to DimArray with respective data.
"""
function loadPreprocData(pathToPreprocDir::String, climateVariables::Vector{String}, diagnostic::String="CLIM", included::Vector{String}=[])
    
    dataAllVars = Dict{String, DimArray}();
    for climVar in climateVariables
        pathToDiagnosticData = joinpath(pathToPreprocDir, climVar * "_" * diagnostic);
        ncFiles = glob("*.nc", pathToDiagnosticData);
        
        data = []
        for file in ncFiles
            addFile = true;
            # if not empty, only includes files that contain names given in 'included'
            if length(included) != 0
                if !all([occursin(name, file) for name in included])
                    addFile = false;
                end
            end
            if addFile
                filename = split(basename(file), ".nc")[1];
                ds = Dataset(file);
                metadata=ds.attrib;

                dsVar = ds[climVar];
                dim1 = Dim{:lon}(collect(dsVar["lon"][:]));
                dim2 = Dim{:lat}(collect(dsVar["lat"][:]));
                dim3 = Dim{:model}([filename]);
                
                mat = Array(dsVar);
                data3d = reshape(mat, size(mat, 1), size(mat, 2), 1)
                push!(data, DimArray(data3d, (dim1, dim2, dim3), metadata=metadata));
            end
        end
        dimData = cat(data..., dims=3)
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
    getModelDistMatrix(data::DimArray)

Computes the area weighted root mean squared error between model predictions for each pair of models. 

data has dimensions: lon, lat, model and contains the data for a single variable.

returns a DimArray of size nxn where n is the number of models in 'data'. 
"""
function getModelDistances(data::DimArray)
    # only take values where none (!) of the models has infinite values!! (Not just the two that are compared to one another)
    nbModels = length(dims(data, :model));
    maskMissing = dropdims(any(ismissing, data, dims=:model), dims=:model);
    latitudes = dims(data, :lat);

    # area weights differ across latitudes (because of the earth being a sphere)
    #areaWeights = cos.(deg2rad.(latitudes))
    #areaWeights = DimArray(areaWeights, Dim{:lat}(Array(latitudes)));
    matrixS = zeros(nbModels, nbModels);
    models = dims(data, :model);
    
    for i = 1 : nbModels
        model_i_name = models[i];
        model_i = data[model=At(model_i_name)];
        model_i[maskMissing .== 1] .= 0;
        
        for j = (i + 1) : nbModels
            model_j_name = models[j];
            model_j = data[model=At(model_j_name)];
            model_j[maskMissing .== 1] .= 0;

            s_ij = areaWeightedMSE(model_i, model_j, maskMissing);
            matrixS[i,j] = s_ij
        end
    end

    symDistMatrix = matrixS .+ matrixS';

    dim = Array(dims(data, :model));
    return DimArray(symDistMatrix, (Dim{:model1}(dim), Dim{:model2}(dim)));
end

"""
    weightDistMatrix(distMatrix::DimArray{T}, weight::T) where T<:Number

Normalizes the matrix by its median and multiplies each entry by 'weight'.
"""
function normalizeAndWeightDistMatrix(distMatrix::Union{DimVector, DimArray{T}}, weight::T=1) where T<:Number
    distMatrix = distMatrix ./ median(distMatrix); 
    weightedDistances = weight .* distMatrix;
    return weightedDistances
end


function normalizeWeightsVariables!(weightsVars::Dict{String, Number})
    total = sum(values(weightsVars))
    for key in keys(weightsVars)
        weightsVars[key] /= total 
    end
end

"""
    calculateWeights(data::DimArray)

Computes an independence weight for each model in 'data'. The weights result from taking the average weighted root mean squared errors between pairs of models for each variable. 

'data' is a Dictionary mapping from the climate variable (e.g. tas) to a DimArray with the corresponding data (lon, lat, model). 
'weightsVars' is a Dictionary mapping from the climate variable to a number which is the weight of how much the respective variable contributes to the computed independenceWeight.

returns a DimArray (model1, model2) with the computed independence weights.
"""
function getIndependenceWeights(data::Dict{String, DimArray}, weightsVars::Dict{String, Number}=nothing)
    if !isnothing(weightsVars)
        normalizeWeightsVariables!(weightsVars);
    end

    weightedDistMatrices = [];
    for climVar in keys(data)
        distances = getModelDistances(data[climVar]);
        weight = ifelse(isnothing(weightsVars), 1, weightsVars[climVar]);
        weightedDistances = normalizeAndWeightDistMatrix(distances, weight);
        push!(weightedDistMatrices, weightedDistances);
    end
    return  reduce(+, weightedDistMatrices);
end

"""
    getModelDataDist(models, observations)

computes the distance (areaweighted rms error) between model predictions and observations. 

"""
function getModelDataDist(models::DimArray, observations::DimArray)      
    distances = [];
    model_names = [];
    for (i, model_i) in enumerate(eachslice(models; dims=:model))
        # indexing using DimensionalData.jl: model_i = models[model=1]
        model_name = dims(models, :model)[i]  # Access the name of the current model

        maskNbMissing = (ismissing.(observations) + ismissing.(model_i)) .> 0; # observations or model is missing (or both)
        maskedObs = deepcopy(observations);
        maskedObs = dropdims(ifelse.(maskNbMissing .> 0, 0, maskedObs), dims=:model);

        maskedModel = ifelse.(maskNbMissing .> 0, 0, model_i);
        mse = areaWeightedMSE(maskedModel, maskedObs, maskNbMissing);

        push!(model_names, model_name);
        push!(distances, mse);
    end
    return DimArray(distances, (Dim{:model}(model_names)))
end


function getPerformanceWeights(modelData::Dict{String, DimArray}, obsData::Dict{String, DimArray}, weightsVars::Dict{String, Number}=nothing)
    if !isnothing(weightsVars)
        normalizeWeightsVariables!(weightsVars);
    end
    variables = keys(modelData);
    weightedDistMatrices = [];

    for climVar in variables
        distances = getModelDataDist(modelData[climVar], obsData[climVar]);
        weight = ifelse(isnothing(weightsVars), 1, weightsVars[climVar]);

        weightedDistances = normalizeAndWeightDistMatrix(distances, weight);
        # just use the actual model name not all the other information stored in the model dimension for averaging over ensemble
        modelNames = map((x) -> split(x, "_")[2], dims(weightedDistances, :model));
        weightedDistances = DimArray(Array(weightedDistances),(Dim{:model}(modelNames)))
        push!(weightedDistMatrices, weightedDistances);
    end
    
    weightedDistances = cat(weightedDistMatrices..., dims=2);
    weightedDistances = DimArray(Array(weightedDistances), (Dim{:model}(Array(dims(weightedDistances, :model))), Dim{:variable}(collect(variables)))); 
    weights = reduce(+, weightedDistances, dims=:variable);

    weightsGrouped = mean.(groupby(weights, dims(weights, :model)));
    indices = .!isnan.(Array(weightsGrouped));

    modelNames = Array(dims(weightsGrouped, :model))
    return DimArray(Array(weightsGrouped)[indices], Dim{:model}(modelNames[indices]))
end
