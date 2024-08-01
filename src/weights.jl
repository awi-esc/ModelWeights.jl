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
        directory = ifelse(isempty(diagnostic), climVar, join([climVar, diagnostic], "_"));
        pathToDiagnosticData = joinpath(pathToPreprocDir, directory);
        ncFiles = glob("*.nc", pathToDiagnosticData);
        
        data = []
        sources = []
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
                ds = Dataset(file);
                meta = ds.attrib;
                if "model_id" in keys(meta)
                    push!(sources, meta["model_id"]);
                else
                    push!(sources, split(filename, "_")[2])
                end
                dsVar = ds[climVar];
                if climVar == "amoc"
                    if "season_number" in keys(ds.dim)
                        dim1 = Dim{:season}(collect(dsVar["season_number"][:]));
                        push!(data, DimArray(Array(dsVar), (dim1), metadata=meta));
                    else
                        push!(data, DimArray(Array(dsVar), (), metadata=meta));
                    end
                else
                    dim1 = Dim{:lon}(collect(dsVar["lon"][:]));
                    dim2 = Dim{:lat}(collect(dsVar["lat"][:]));          
                    push!(data, DimArray(Array(dsVar), (dim1, dim2), metadata=meta));
                end
            end
        end
        dimData = cat(data...; dims = (Dim{:model}(sources)));
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
    getIndependenceWeights(data::Dict{String, DimArray}, weightsVars::Dict{String, Number}=nothing)

Computes an independence weight for each model in 'data'. The weights result from taking the average weighted root mean squared errors between pairs of models for each variable. 

'data' is a Dictionary mapping from the climate variable (e.g. tas) to a DimArray with the corresponding data (lon, lat, model). 
'weightsVars' is a Dictionary mapping from the climate variable to a number which is the weight of how much the respective variable contributes to the computed independenceWeight.

returns a DimArray (model1, model2, variable) with the computed independence weights, seperately for each considered variable.
"""
function getIndependenceWeights(data::Dict{String, DimArray}, weightsVars::Dict{String, Number}=nothing)
    weights = copy(weightsVars);
    if !isnothing(weightsVars)
        normalizeWeightsVariables!(weights);
    end
    weightedDistMatrices = [];
    variables = keys(data)
    for climVar in variables
        distances = getModelDistances(data[climVar]);
        weight = ifelse(isnothing(weights), 1, weights[climVar]);
        weightedDistances = normalizeAndWeightDistMatrix(distances, weight);
        push!(weightedDistMatrices, weightedDistances);
    end
    weightsByVars = cat(weightedDistMatrices..., dims = Dim{:variable}(collect(variables)));
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
    getPerformanceWeights(modelData::Dict{String, DimArray}, obsData::Dict{String, DimArray}, weightsVars::Dict{String, Number}=nothing)

Computes a performance weight for each model in 'modelData'. The weights result from taking the average weighted root mean squared errors between each model and the observational data.
'mdoelData' and 'obsData' are Dictionaries mapping from the climate variable (e.g. tas) to a DimArray with the corresponding data (lon, lat, model). For the observational data, 
the third dimension (model) just contains a single entry (e.g. ERA5).
'weightsVars' is a Dictionary mapping from the climate variable to a number which is the weight of how much the respective variable contributes to the computed performanceWeight.

returns a DimArray (model1) with the computed independence weights, seperately for each considered variable.
"""
function getPerformanceWeights(modelData::Dict{String, DimArray}, 
                               obsData::Dict{String, DimArray},
                               weightsVars::Dict{String, Number}=nothing
                               )
    weights = copy(weightsVars);
    if !isnothing(weights)
        normalizeWeightsVariables!(weights);
    end
    variables = keys(modelData);
    weightedDistMatrices = [];

    for climVar in variables
        distances = getModelDataDist(modelData[climVar], obsData[climVar]);
        weight = ifelse(isnothing(weights), 1, weights[climVar]);

        weightedDistances = normalizeAndWeightDistMatrix(distances, weight);
        # just use the actual model name not all the other information stored in the model dimension for averaging over ensemble
        #modelNames = map((x) -> split(x, "_")[2], dims(weightedDistances, :model));
        #weightedDistances = DimArray(Array(weightedDistances),(Dim{:model}(modelNames)))
        push!(weightedDistMatrices, weightedDistances);
    end
    
    # TODO: here metadata has to be combined correctly
    weightsByVar = cat(weightedDistMatrices..., dims = Dim{:variable}(collect(variables)));
    return weightsByVar
end

function summarizeWeightsAcrossVars(weightsByVar::DimArray)
    summarizedWeights = reduce(+, weightsByVar, dims=:variable);
    return dropdims(summarizedWeights, dims=:variable)
end


function averageEnsembleVector(weights::DimArray)
    # TODO: here metadata gets lost, just groups are retained
    weightsGrouped = mean.(groupby(weights, dims(weights, :model)));
    indices = .!isnan.(Array(weightsGrouped));
    modelNames = Array(dims(weightsGrouped, :model))
    # TODO: add metadata!
    weightsByEnsemble =  DimArray(Array(weightsGrouped)[indices], Dim{:model}(modelNames[indices]));
    return weightsByEnsemble
end



function averageEnsembleMatrix(distances::DimArray)
    # TODO: add metadata
    models = dims(distances, :model1);
    modelsUnique = unique(models);
    
    ensemble = DimArray(zeros(length(modelsUnique), length(modelsUnique)), (Dim{:model1}(modelsUnique), Dim{:model2}(modelsUnique)));
    for model in modelsUnique
        avg = mean(distances[models.!=model, models.==model], dims=:model2)
        dimModel2 = [model];  
        avg = set(avg, :model2 => dimModel2)
        dimModel1 = Array(dims(avg, :model1))
        ensemble[model1=At(dimModel1), model2=At(dimModel2)] = avg
    end
    return ensemble
end

"""
    combineWeights(performanceWeights::DimArray, independenceWeights::DimArray, sigmaD::Number=0.5, sigmaS::Number=0.5)

Combines the RMSEs between pairs of models and the RMSEs between each model and the data into a set of normalized weights, 
one for each model. 

Note that, the given weights are summarized wrt ensemble models, s.t. for each ensemble there is only one value. 

'sigmaD': Nb. btw. 0 and 1; contribution of performance metric 
'sigmaS': Nb. btw. 0 and 1; contribution of independence metric
"""
function combineWeights(performanceWeights::DimArray, independenceWeights::DimArray, sigmaD::Number, sigmaS::Number)
    
    wP = averageEnsembleVector(performanceWeights);
    wI = averageEnsembleMatrix(independenceWeights);
    performance = exp.(-(wP ./ sigmaD).^2);
    # note: (+1 in eq. in paper is from when model is compared to itself since exp(0)=1)
    independence = dropdims(sum(exp.(-(wI ./ sigmaS).^2), dims=:model2), dims=:model2)
    independence = DimArray(Array(independence), (Dim{:model}(Array(dims(independence, :model1)))))

    weights = performance ./ independence;
    weightsNormalized = weights ./ sum(weights);
    return weightsNormalized
end

"""
    getWeights(modelData::DimArray, obsData::DimArray, sigmaD::Number=0.5, sigmaS::Number=0.5, 
               weightsPerform::Dict{String, Number}=nothing, 
               weightsIndep::Dict{String, Number}=nothing
               )

    Computes weight for each model in multi-model ensemble according to approach from 
    Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner, Anna L. Merrifield, Ruth Lorenz, and Reto Knutti.
    “Reduced Global Warming from CMIP6 Projections When Weighting Models by Performance and Independence.” 
    Earth System Dynamics 11, no. 4 (November 13, 2020): 995–1012. https://doi.org/10.5194/esd-11-995-2020.

"""
function getWeights(modelData::Dict{String, DimArray}, 
                    obsData::Dict{String, DimArray}, 
                    sigmaD::Number=0.5, 
                    sigmaS::Number=0.5,
                    weightsPerform::Dict{String, Number}=nothing, 
                    weightsIndep::Dict{String, Number}=nothing
    )
    wP = getPerformanceWeights(modelData, obsData, weightsPerform);
    wI = getIndependenceWeights(modelData, weightsIndep);
    weights = combineWeights(summarizeWeightsAcrossVars(wP), summarizeWeightsAcrossVars(wI), sigmaD, sigmaS)
    return weights
end