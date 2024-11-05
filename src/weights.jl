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
    generalizedDistances(
        modelData::Dict{String, DimArray},
        dist_type::String;
        weightsVars::Dict{String, Number}=Dict{String, Number}(),
        obsData::Dict{String, DimArray}
    )

Compute the generalized distance for each model and climate variable in 'modelData' with respect to another model (dist_type="independence") or
to observational data (dist_type = "performance").

# Arguments
- `modelData`: keys are climate variables, values are DimArrays with dimensions 'lon', 'lat' and 'model' (dist_type=performance) or 'model1' and 'model2' (dist_type=independence) 
- `obsData`:  keys are climate variables, values are DimArrays with dimensions 'lon', 'lat', 'model'
- `weightsVars`: keys are climate variables, values are the weight of how much the respective variable contributes to the computed weight

# Returns 
- `weightsByVar::DimArray`: performance or independence weights for each considered variable
"""
function generalizedDistances(
    modelData::Dict{String, DimArray}, 
    dist_type::String;
    weightsVars::Dict{String, Number}=Dict{String, Number}(),
    obsData::Dict{String, DimArray}=Dict{String, DimArray}()
)
    variables = collect(keys(weightsVars))
    weights = copy(weightsVars);
    if !isempty(weights)
        normalizeWeightsVariables!(weights);
    else
        weights = Dict(zip(variables, ones(length(variables))));
    end
    weightedDistMatrices = [];

    meta = Dict{String, Union{String, Array, Dict}}();
    for climVar in variables
        modelDict = filter(((k,v),) -> startswith(k, climVar * "_"), modelData)
        obsDataDict = filter(((k,v),) -> startswith(k, climVar * "_"), obsData)
        if length(modelDict) > 1
            @warn "more than one dataset for computing generalizedDistances for climate variable $climVar in model data, first is taken!"
        end
        model = first(values(modelDict))
        
        if dist_type == "performance"
            if length(obsDataDict) > 1
                @warn "more than one dataset for computing generalizedDistances for climate variable $climVar in observational data, first is taken!"
            end
            obs = first(values(obsDataDict))
            distances = getModelDataDist(model, obs)
        elseif dist_type == "independence"
            distances = getModelDistances(model)
        else
            throw(ArgumentError("Argument 'dist_type' in generalizedDistances must be one of: 'performance', 'independence'."))
        end
        weightedNormalizedDistances = normalizeAndWeightDistances(distances, weights[climVar]);
        push!(weightedDistMatrices, weightedNormalizedDistances);

        metadata = deepcopy(model.metadata);
        meta_shared = getMetadataSharedAcrossModelsAndModelNames(metadata);
        meta_shared["variables"] = climVar
        meta = mergewith(appendValuesDicts, meta, meta_shared);
    end
    weightsByVar = cat(weightedDistMatrices..., dims = Dim{:variable}(collect(variables)));
    weightsByVar = rebuild(weightsByVar; metadata = meta);
    return weightsByVar
end


""" 
    averageEnsembleVector(data::DimArray, updateMeta::Bool)

For each model and variable (if several given), compute the mean across all
ensemble members of that model.

# Arguments:
- `data`: a DimArray with dimension 'model' and possibly 'variable'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
If true attribute ensemble_indices_map is set in metadata.
Set to false if vectors refer to different variables for instance. 

# Return:
- DimArray with dimensions 'model' and possibly 'variable'
"""
function averageEnsembleVector(data::DimArray, updateMeta::Bool)
    grouped = nothing;
    data = renameModelDimsFromMemberToEnsemble(data, ["model"])
    
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
If true attribute ensemble_indices_map is set in metadata. 
Set to false if vectors refer to different variables for instance. 
"""
function averageEnsembleMatrix(data::DimArray, updateMeta::Bool)
    data = renameModelDimsFromMemberToEnsemble(data, ["model1", "model2"])
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
    performanceParts(generalizedDistances::DimArray, sigmaD::Number)

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
    model_data::Data, 
    obs_data::Data,
    config_weights::ConfigWeights
)

    obs_keys = map(x -> (var = x.variable, stat = x.statistic), obs_data.ids)
    model_keys = map(x -> (var = x.variable, stat = x.statistic), model_data.ids)
    if obs_keys != model_keys
        msg = "Variable+Diagnostic combinations are not the same for observational and model data! Obs: $obs_keys , Model: $model_keys"
        throw(ArgumentError(msg))
    end

    diagnostics = unique(map(x -> x.stat, obs_keys))
    Di_all = []
    Sij_all = []
    Di = []
    Sij = []
    for diagnostic in diagnostics
        models = filter(((k, v),) -> occursin("_" * diagnostic * "_", k), model_data.data)
        obsData = filter(((k, v),) -> occursin("_" * diagnostic * "_", k), obs_data.data)

        weights_perform = filter(((k, v),) -> occursin("_" * diagnostic, k), config_weights.performance)
        weights_indep = filter(((k, v),) -> occursin("_" * diagnostic, k), config_weights.independence)
        # weights dictionary for this particular diagnostic shall map just from variable to value
        wP = Dict{String, Number}()
        wI = Dict{String, Number}()
        for key in keys(weights_perform)
            wP[split(key, "_")[1]] = weights_perform[key]
            wI[split(key, "_")[1]] = weights_indep[key]
        end
        generalizedDistsPerform  = generalizedDistances(
            models, "performance"; weightsVars=wP, obsData=obsData
        )
        push!(Di_all, generalizedDistsPerform)
        push!(Di, reduceGeneralizedDistancesVars(generalizedDistsPerform))

        generalizedDistsIndep = generalizedDistances(models, "independence"; weightsVars=wI);
        push!(Sij_all, generalizedDistsIndep)
        push!(Sij, reduceGeneralizedDistancesVars(generalizedDistsIndep))
    end
    # put into DimArray
    Di_all = cat(Di_all..., dims = Dim{:diagnostic}(collect(diagnostics)));
    Sij_all = cat(Sij_all..., dims = Dim{:diagnostic}(collect(diagnostics)));
    
    # summarized across variables (done in for loop) and diagnostics
    Di = cat(Di..., dims = Dim{:diagnostic}(collect(diagnostics)));
    Sij = cat(Sij..., dims = Dim{:diagnostic}(collect(diagnostics)));
    Di = dropdims(reduce(+, Di, dims=:diagnostic), dims=:diagnostic)
    Sij = dropdims(reduce(+, Sij, dims=:diagnostic), dims=:diagnostic)

    return (Di_all, Sij_all, Di, Sij)
end


"""
    computeWeightedAvg(data::DimArray, w::Union{DimArray, Nothing}=nothing)

Compute the average values for each (lon,lat) grid point in 'data_var', weighted
by weights 'w'. Members of the same ensemble are averaged before. If no weight
vector is provided, unweighted average is computed.

# Arguments:
- `data`: DimArray with dimensions lon, lat, model
- `w`: DimArray with dimension 'model'
"""
function computeWeightedAvg(
    data::DimArray; weights::Union{DimArray, Nothing}=nothing
)
    data = deepcopy(data)
    if isnothing(weights)
        # make sure that the number of ensemble members per model is considered
        ensemble_names = data.metadata["ensemble_names"]
        indices = []
        for model in unique(ensemble_names)
            push!(indices, findall(x -> x==model, ensemble_names))
        end
        n_ensembles = length(indices)
        w = []
        for positions in indices
            n_members = length(positions)
            for _ in range(1, n_members)
                push!(w, (1/n_ensembles) *  (1/n_members))
            end
        end
        weights = DimArray(w, Dim{:model}(Array(dims(data, :model))))
    else
        # weights may have been computed wrt a different set of variables as we use here, 
        # so the list of models for which weights have been computed may be shorter 
        # than the models of the given data (for the same reference period).
        if sort(dims(weights, :model)) != sort(models)
            @warn "Mismatch between models that weights were computed for and models in the data."
        end
        data = data[model = Where(m -> m in dims(weights, :model))]
        n_models_data = length(dims(data, :model))
        n_weights = length(weights)
        if n_models_data != n_weights
            msg = "nb of models for observational and model predictions does not match: ";
            msg2 = "weights: " * string(n_weights) * " , data: " * string(n_models_data);
            throw(ArgumentError(msg * msg2))
        end
    end
    @assert isapprox(sum(weights), 1; atol=10^-4)
    for m in dims(data, :model)
        data[model = At(m)] = data[model = At(m)] .* weights[model = At(m)]
    end

    weighted_avg = dropdims(reduce(+, data, dims=:model), dims=:model)
    return weighted_avg
end


function computeWeights(
    model_data::Data, obs_data::Data, config_weights::ConfigWeights
)
    Di_all, Sij_all, Di, Sij = overallGeneralizedDistances(
        model_data, obs_data, config_weights
    );
    performances = performanceParts(Di, config_weights.sigma_performance);
    independences = independenceParts(Sij, config_weights.sigma_independence);
    weights = performances ./ independences;
    weights = weights ./ sum(weights);
    weights.metadata["name_ref_period"] = config_weights.ref_period
    
    wP = performances ./ sum(performances)
    wI = independences ./ sum(independences)
    #w = wP./wI # just for sanity check

    return ClimwipWeights(
        performance_all = Di_all, 
        independence_all = Sij_all, 
        Di = Di,
        Sij = Sij,
        wP = wP,
        wI = wI,
        w =  weights
        #overall = w./sum(w), # just for sanity check
    )
end


function makeWeightPerEnsembleMember(weights::DimArray)
    full_model_names = weights.metadata["full_model_names"]
    models = weights.metadata[getCMIPModelsKey(weights.metadata)]
    nb_ensemble_members = [];
    for model in models
        n = count(member -> startswith(member, model), full_model_names)
        push!(nb_ensemble_members, n)
    end
    weights_all_members = []
    for (idx, w) in enumerate(weights)
        n_members = nb_ensemble_members[idx]
        for _ in range(1, n_members)
            push!(weights_all_members, w/nb_ensemble_members[idx])
        end
    end
    return DimArray(weights_all_members, Dim{:model}(full_model_names), metadata = weights.metadata)
end


function logWeights(metadata_weights)
    models = metadata_weights["full_model_names"];
    @info "Nb included models (without ensemble members): " length(metadata_weights["ensemble_names"])
    foreach(m -> @info(m), models)
end


function applyWeights(data::Data, weights::ClimwipWeights)
    
end