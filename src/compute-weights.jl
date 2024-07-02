using NCDatasets
using DimensionalData
using Glob
using Statistics

"""
    loadPreprocData(pathToPreprocData, climateVar[, diagnostic=CLIM, dataIncluded=[])

Loads all .nc files (each model is a different .nc file) for variable climateVar and 'diagnostic' inside the respective subdirectory (climateVar_diagnostic)
of the directory located at 'path2PreprocData' into a single DimArray with dimensions lon (longitude), lat (latitude) and model. 

If length(dataIncluded) != 0 only those .nc files are considered whose filenames contain all elements of dataIncluded, 
e.g. if dataIncluded=['ERA5'] only ERA5 data will be included (files with ERA5 in their filename).

returns a dictionary from climateVariable to DimArray with respective data.
"""
function loadPreprocData(pathToPreprocDir::String, climateVariables::Vector{String}, diagnostic::String, included::Vector{String}=[])
    
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
    areaWeights = cos.(deg2rad.(latitudes))
    dim = Dim{:lat}(Array(latitudes));
    areaWeights = DimArray(areaWeights, dim);
    matrixS = zeros(nbModels, nbModels);
    # apply mask: ignore entries where at least one model has infinite value
    models = dims(data, :model);
    
    for i = 1 : nbModels
        model_i_name = models[i];
        model_i = data[model=At(model_i_name)];
        model_i[maskMissing .== 1] .= 0;
        
        for j = (i + 1) : nbModels
            model_j_name = models[j];
            model_j = data[model=At(model_j_name)];
            model_j[maskMissing .== 1] .= 0;

            sqDiff = (model_i .- model_j).^2;
            values = areaWeights' .* sqDiff;
            s_ij = sqrt(sum(values));
            matrixS[i,j] = s_ij
        end
    end

    symDistMatrix = matrixS .+ matrixS';

    dim = Array(dims(data, :model));
    return DimArray(symDistMatrix, (Dim{:model1}(dim), Dim{:model2}(dim)));
end

"""
    weightDistMatrix(distMatrix::DimArray{T}, weight::T) where T<:Number

Normalizes the matrix by its mean and multiplies each entry by 'weight'.
"""
function normalizeAndWeightDistMatrix(distMatrix::DimArray{T}, weight::T=1) where T<:Number
    distMatrix = distMatrix ./ median(distMatrix);     # normalize distance matrix
    weightedDistances = weight .* distMatrix;
    return weightedDistances
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
        # normalize weights for contribution of variables to make them sum up to 1
        total = sum(values(weightsVars))
        for key in keys(weightsVars)
            weightsVars[key] /= total 
        end
    end
    weightedDistMatrices = [];
    for climVar in keys(data)
        distances = getModelDistances(data[climVar]);
        if isnothing(weightsVars)
            weight = 1;
        else
            weight = weightsVars[climVar];
        end
        weightedDistances = normalizeAndWeightDistMatrix(distances, weight);
        push!(weightedDistMatrices, weightedDistances);
    end
    independenceWeights = reduce(+, weightedDistMatrices);
    return independenceWeights
end

