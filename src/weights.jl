using NCDatasets
using DimensionalData
using Statistics
using LinearAlgebra

"""
    areaWeightedRMSE(m1::DimArray, m2::DimArray, mask::DimArray)

Compute the area weighted (cosine of latitudes in radians) root mean squared 
error between two matrices. 

# Arguments:
- `m1`: has dimensions 'lon', 'lat'
- `m2`: has dimensions 'lon', 'lat'
- `mask`: has values 0,1. Locations where mask is 1 are ignored, i.e. they get a weight of 0!

# Return:
- single value, area-weighted root mean squared error
"""
function areaWeightedRMSE(m1::DimArray, m2::DimArray, mask::DimArray)
    latitudes = dims(m1, :lat);
    areaWeights = cos.(deg2rad.(latitudes));
    areaWeights = DimArray(areaWeights, Dim{:lat}(Array(latitudes)));

    sqDiff = (m1 .- m2).^2;
    weightedValues = areaWeights' .* sqDiff;

    weightMatrix = repeat(areaWeights', length(DimensionalData.dims(m1, :lon)), 1);  
    weightMatrix = ifelse.(mask .== 1, 0, weightMatrix); 
    normalization = sum(weightMatrix);

    return sqrt(sum(skipmissing(weightedValues)./normalization));
end

"""
    getModelDistMatrix(modelData::DimArray)

Computes the area weighted root mean squared error between model predictions for each pair of models. 

# Arguments:
- `modelData` has dimensions 'lon', 'lat', 'model' and contains the data for a single climate variable.

# Return:
- Symmetrical matrix (::DimArray) of size nxn where n is the number of models in 'modelData'. 
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
            s_ij = areaWeightedRMSE(model_i, model_j, maskMissing);
            matrixS[i, idx] = s_ij
        end
    end
    # make symmetrical matrix
    symDistMatrix = matrixS .+ matrixS';
    dim = Array(dims(modelData, :model));
    return DimArray(symDistMatrix, (Dim{:model1}(dim), Dim{:model2}(dim)))
end


"""
    normalizeAndWeightDistances(distMatrix::Union{DimVector, DimArray}, weight::T=1) where T<:Real

Normalize 'distMatrix' by its median and multiplies each entry by 'weight'.
"""
function normalizeAndWeightDistances(distMatrix::Union{DimVector, DimArray}, weight::T=1) where T<:Real
    distMatrix = distMatrix ./ median(distMatrix); 
    return weight .* distMatrix;
end


"""
    normalizeWeightsVariables!(weightsVars::Dict{String, Number})

Modify weights for each variable by normalizing them s.t. they sum up to 1.
"""
function normalizeWeightsVariables!(weights::Dict{String, Number})
    total = sum(values(weights))
    for key in keys(weights)
        weights[key] /= total 
    end
end


"""
    generalizedDistancesIndependence(
        data::Dict{String, DimArray}, 
        weightsVars::Dict{String, Number}=Dict{String, Number}()
    )

Compute an independence weight for each model and variable in 'data'. The 
weights are the average weighted root mean squared errors between pairs of 
models (for each variable). 

# Arguments:
- `data`: keys are climate variables (e.g. tas), values are DimArrays with the
corresponding data and dimensions 'lon', 'lat', 'model'. 
- `weightsVars`: keys are climate variables, values are the weights of how much
the respective variable contributes to the computed independenceWeight.

# Return:
- `weightsByVars`: a DimArray with dimensions 'model1', 'model2', 'variable'
with the computed independence weights, seperately for each considered variable.
"""
function generalizedDistancesIndependence(
    data::Dict{String, DimArray}, 
    weightsVars::Dict{String, Number}=Dict{String, Number}()
)
    variables = keys(data);
    weights = copy(weightsVars);
    if !isempty(weights)
        normalizeWeightsVariables!(weights);
    else
        weights = Dict(zip(variables, ones(length(variables))));
    end

    weightedDistMatrices = [];
    meta = Dict{String, Array}();
    # Models are/must be the same across variables! 
    for climVar in variables
        distances = getModelDistances(data[climVar]);
        weight = ifelse(isempty(weights), 1, weights[climVar]);        
        weightedNormalizedDistances = normalizeAndWeightDistances(distances, weight);
        push!(weightedDistMatrices, weightedNormalizedDistances);
        
        metadata = deepcopy(data[climVar].metadata);
        meta_shared = getSharedMetadataAndModelNames(metadata);
        meta_shared["variables"] = climVar
        meta = mergewith(appendValuesDicts, meta, meta_shared);
    end
    weightsByVars = cat(weightedDistMatrices..., dims = Dim{:variable}(collect(variables)));
    meta["indices_map"] = getIndicesMapping(collect(variables))
    weightsByVars = rebuild(weightsByVars; metadata = meta);
    return weightsByVars
end

"""
    getModelDataDist(models::DimArray, observations::DimArray)

Compute the distance (area-weighted rmse) between model predictions and observations. 

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
        mse = areaWeightedRMSE(maskedModel, maskedObs, maskNbMissing);

        push!(model_names, model_name);
        push!(distances, mse);
    end
    return DimArray(distances, (Dim{:model}(model_names)), metadata=models.metadata)
end


"""
    reduceGeneralizedDistancesVars(normalizedWeightedDistsByVar::DimArray)

Compute the generalized distance as the sum across normalized and weighted
values for each variable (and diagnostic). 
"""
function reduceGeneralizedDistancesVars(normalizedWeightedDistsByVar::DimArray)
    generalizedDist = reduce(+, normalizedWeightedDistsByVar, dims=:variable)
    return dropdims(generalizedDist, dims=:variable)
end


"""
    generalizedDistancesPerformance(
        modelData::Dict{String, DimArray}, 
        obsData::Dict{String, DimArray}, 
        weightsVars::Dict{String, Number}=Dict{String, Number}()
    )

Compute a performance weight for each model and climate variable in 'modelData'.
The weights are the average weighted root mean squared errors between each model and the observational data.

# Arguments
- `modelData`: keys are climate variables, values are DimArrays with dimensions 'lon', 'lat', 'model' 
- `obsData`:  keys are climate variables, values are DimArrays with dimensions 'lon', 'lat', 'model'
- `weightsVars`: keys are climate variables, values are the weight of how much the respective variable contributes to the computed performanceWeight

# Returns 
- `weightsByVar::DimArray`: performance weights for each considered variable
"""
function generalizedDistancesPerformance(
    modelData::Dict{String, DimArray}, 
    obsData::Dict{String, DimArray},
    weightsVars::Dict{String, Number}=Dict{String, Number}()
)
    variables = keys(modelData);
    weights = copy(weightsVars);
    if !isempty(weights)
        SimilarityWeights.normalizeWeightsVariables!(weights);
    else
        weights = Dict(zip(variables, ones(length(variables))));
    end
    weightedDistMatrices = [];

    meta = Dict{String, Union{String, Array, Dict}}();
    for climVar in variables
        distances = SimilarityWeights.getModelDataDist(modelData[climVar], obsData[climVar]);
        weight = ifelse(isempty(weights), 1, weights[climVar]);
        weightedNormalizedDistances = SimilarityWeights.normalizeAndWeightDistances(distances, weight);
        push!(weightedDistMatrices, weightedNormalizedDistances);

        metadata = deepcopy(modelData[climVar].metadata);
        meta_shared = SimilarityWeights.getSharedMetadataAndModelNames(metadata);
        meta_shared["variables"] = climVar
        meta = mergewith(SimilarityWeights.appendValuesDicts, meta, meta_shared);
    end
    weightsByVar = cat(weightedDistMatrices..., dims = Dim{:variable}(collect(variables)));
    meta["indices_map"] = getIndicesMapping(collect(variables))
    weightsByVar = rebuild(weightsByVar; metadata = meta);
    return weightsByVar
end


""" 
    averageEnsembleVector!(data::DimArray, updateMeta::Bool)

For each model and variable (if several given), compute the mean across all
ensemble members of that model.

# Arguments:
- `data`: a DimArray with dimension 'model' and possibly 'variable'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
If true attribute indices_map is set in metadata. Set to false if vectors refer to 
different variables for instance. 

# Return:
- DimArray with dimensions 'model' and possibly 'variable'
"""
function averageEnsembleVector(data::DimArray, updateMeta::Bool)
    grouped = nothing;
    if !hasdim(data, :variable)
        grouped = groupby(data, :model=>identity);
        averages = map(d -> mean(d, dims=:model), grouped);
        models = Array(dims(averages, :model));
        combined = cat(averages..., dims=(Dim{:model}(models)));
    else
        variables = Array(dims(data, :variable));
        results = [];
        for climVar in variables
            grouped = groupby(data[variable = Where(x -> x == climVar)], :model=>identity);
            averages = map(d->mean(d, dims=:model), grouped);
            models = Array(dims(averages, :model));
            avg = cat(averages..., dims=(Dim{:model}(models)));
            push!(results, avg);
        end
        models = unique(Array(dims(data, :model)));
        combined = cat(results...; dims=(Dim{:variable}(variables)));
    end

    if updateMeta
        meta = updateGroupedDataMetadata(data.metadata, grouped)
    else 
        meta = data.metadata
    end
    combined = rebuild(combined; metadata = meta);
    return combined;
end


"""
    averageEnsembleMatrix(data::DimArray, updateMeta::Bool)

Compute the average weights across all members of each model ensemble for each
given variable.

# Arguments:
- `data`: DimArray with dimensions 'model1', 'model2' and possibly 'variable'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
If true attribute indices_map is set in metadata. Set to false if vectors refer to 
different variables for instance. 
"""
function averageEnsembleMatrix(data::DimArray, updateMeta::Bool)
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
    if updateMeta
        meta = updateGroupedDataMetadata(data.metadata, grouped)
    else 
        meta = data.metadata
    end
    combined = rebuild(combined; metadata = meta);
    return combined
end


function averageEnsembleMembers(generalizedDistances::DimArray, updateMeta::Bool)
    if hasdim(generalizedDistances, Symbol("model1"))
        distances = averageEnsembleMatrix(generalizedDistances, updateMeta);
    else 
        distances = averageEnsembleVector(generalizedDistances, updateMeta);
    end
    return distances
end


"""
    performanceParts(generalizedDistances::DimArray)

Compute the performance part (numerator) of the overall weight for each model.

# Arguments:
- `generalizedDistances`: has a single dimension 'model' and contains D_i for 
each model i
- `sigmaD`: free model parameter for impact of performance weights
"""
function performanceParts(generalizedDistances::DimArray, sigmaD::Number)
    distances = averageEnsembleMembers(generalizedDistances, false);
    return exp.(-(distances ./ sigmaD).^2);
end


"""
    independenceParts(generalizedDistances::DimArray, sigmaS::Number)

Compute the independence part (denominator) of the overall weight for each model.

# Arguments:
- `generalizedDistances`: has two dimensions, 'model1' and 'model2' and
contains S_{i,j} for each model pair 
- `sigmaS`: free model parameter for impact of independence weights
"""
function independenceParts(generalizedDistances::DimArray, sigmaS::Number)
    distances = averageEnsembleMembers(generalizedDistances, false);
    # note: (+1 in eq. in paper is from when model is compared to itself since exp(0)=1)
    indep_parts = sum(exp.(-(distances ./ sigmaS).^2), dims=:model2);
    indep_parts = dropdims(indep_parts, dims=:model2)
    indep_parts = set(indep_parts, :model1 => :model);
    return indep_parts
end


function overallGeneralizedDistances(
    modelData, 
    obsData, 
    weightsPerform, 
    weightsIndep
)
    Di_all_diagnostics = []
    Sij_all_diagnostics = []
    diagnostics = keys(modelData)
    for diagnostic in diagnostics
        generalizedDistsPerform  = generalizedDistancesPerformance(
            modelData[diagnostic], obsData[diagnostic], weightsPerform
        );
        Di = reduceGeneralizedDistancesVars(generalizedDistsPerform);
        push!(Di_all_diagnostics, Di)

        generalizedDistsIndep = generalizedDistancesIndependence(
            modelData[diagnostic], weightsIndep
        );
        Sij =  reduceGeneralizedDistancesVars(generalizedDistsIndep);
        push!(Sij_all_diagnostics, Sij)
    end
    Di = cat(Di_all_diagnostics..., dims = Dim{:diagnostic}(collect(diagnostics)));
    Sij = cat(Sij_all_diagnostics..., dims = Dim{:diagnostic}(collect(diagnostics)));
    
    Di = dropdims(reduce(+, Di, dims=:diagnostic), dims=:diagnostic)
    Sij = dropdims(reduce(+, Sij, dims=:diagnostic), dims=:diagnostic)
    return (Di, Sij)
end


"""
    computeWeightedAvg(data_var::DimArray, w::Union{DimArray, Nothing}=nothing)

Compute the average values for each (lon,lat) grid point in 'data_var', weighted
by weights 'w'. Members of the same ensemble are averaged before. If no weight
vector is provided, unweighted average is computed.

# Arguments:
- `data_var`: DimArray with dimensions lon, lat, model
- `w`: DimArray with dimension 'model'
"""
function computeWeightedAvg(data_var::DimArray, w::Union{DimArray, Nothing}=nothing)
    data = averageEnsembleVector(data_var, true);
    models = dims(data, :model);
    
    if isnothing(w)
        w = DimArray(ones(length(models))./length(models), Dim{:model}(Array(models)))
    else
        # weights may have been computed wrt a different set of variables as we use here, 
        # so the list of models for which weights have been computed may be shorter 
        # than the models of the given data (for the same reference period).
        data = data[model = Where(m -> m in dims(w, :model))]
        models = dims(data, :model)
        if length(models) != length(w)
            msg = "nb of models for observational and model predictions does not match: ";
            msg2 = "weights: " * length(w) * " , data: " * length(models);
            throw(ArgumentError(msg * msg2))
        end
    end

    for m in models
        data[model = At(m)] = data[model = At(m)] .* w[model = At(m)]
    end

    weighted_avg = dropdims(reduce(+, data, dims=:model), dims=:model)
    return weighted_avg
end


function computeWeights(
    modelData, 
    obsData, 
    weights_vars_perform, 
    weights_vars_indep, 
    sigmaD, 
    sigmaS
)
    Di, Sij = overallGeneralizedDistances(
        modelData, obsData, weights_vars_perform, weights_vars_indep
    );
    performances = performanceParts(Di, sigmaD);
    independences = independenceParts(Sij, sigmaS);
    weights = performances ./ independences;
    weights = weights ./ sum(weights);
    return weights
end

function logWeights(weights_metadata)
    models = weights_metadata["full_model_names"];
    model_key = getCMIPModelsKey(weights_metadata);
    @info "Nb included models (without ensemble members): " length(weights_metadata[model_key])
    foreach(m -> @info(m), models)
end