using NetCDF
using Statistics
using DataStructures

include(joinpath(@__DIR__, "..", "src", "load-data-utils.jl"));

# 1. Preprocessed data for weighted_temperature_graph diagnostics
path2preprocWeightedTempGraph = joinpath(PATH_TO_PREPROC_DIR, "weighted_temperature_graph", "tas")
preprocTempGraph, modelFileNames = loadNCdataInDir(path2preprocWeightedTempGraph, "tas", []);

# get the stored data from work-directory that was used in diagnostics
path2TempAnomalies = joinpath(PATH_TO_WORK_DIR,  "weighted_temperature_graph", "weighted_temperature_graph", "temperature_anomalies.nc");
tempAnomalies = NetCDF.ncread(path2TempAnomalies, "tas")

# the data in the preproc directory should be identical to the data stored in work directory
@assert tempAnomalies == preprocTempGraph
# then, the unweighted averages should also be identical
unweightedAvg = Statistics.mean(tempAnomalies, dims = 2);
unweightedAvgFromPreproc = Statistics.mean(preprocTempGraph, dims = 2);
@assert unweightedAvgFromPreproc == unweightedAvg


# 2. calculate_weights_climwip
# 2.1. Independence: Model-Model Distance matrices 
PATH_TO_PREPROC_WEIGHTS_DIR = joinpath(PATH_TO_PREPROC_DIR, "calculate_weights_climwip");

#################### do some checks ###########################################
# look at output from single model
pathSingleModel = joinpath(PATH_TO_PREPROC_DIR, "calculate_weights_climwip", "tas_CLIM", "CMIP5_ACCESS1-3_Amon_historical-rcp85_r1i1p1_tas_1995-2014.nc");
pathSingleModel = joinpath(PATH_TO_PREPROC_DIR, "calculate_weights_climwip", "tas_CLIM", "CMIP5_BNU-ESM_Amon_historical-rcp85_r1i1p1_tas_1995-2014.nc");
dataSingleModel = NetCDF.ncread(pathSingleModel, "tas");
size(dataSingleModel)
timeSingleModel = NetCDF.ncread(pathSingleModel, "time");
timeBndsSingleModel = NetCDF.ncread(pathSingleModel, "time_bnds");

@assert timeSingleModel[1] == Statistics.mean(timeBndsSingleModel)

# observational data
# mean across time bounds had already been computed, therefore for each lat x lon - pair just one value; result is matrix of size (lon x lat) 
preprocERA5Tas, _ = loadNCdataInDir(joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, "tas_CLIM"), "tas", ["ERA5"]);
size(preprocERA5Tas)
# ncols = size(dataSingleModel)[2]
# nbModels = Int(size(preprocTas)[2] / ncols)
latitudesAll, _ = loadNCdataInDir(joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, "tas_CLIM"), "lat", ["CMIP"])
latitudesAll, _ = loadNCdataInDir(joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, "pr_CLIM"), "lat", ["CMIP"])
## latitudes should be identical across models
latitudes = NetCDF.ncread(pathSingleModel, "lat")
###############################################################

function loadModelData(climateVar::String)
    climateVarDir = climateVar * "_CLIM";
    preprocData, filenames = loadNCdataInDir(joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, climateVarDir), climateVar, ["CMIP"], false);
    latitudes = NetCDF.ncread(joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, climateVarDir, filenames[1]), "lat");
    preprocData3d = cat(preprocData..., dims=3);
    return preprocData3d, latitudes
end

function getInfMask(data3d, dim=3)
    n = size(data3d)[3];
    maskNotInf = data3d .!= 1.0f20;
    maskNoneInf = sum(maskNotInf, dims=dim) .== n;  
    maskNoneInf = dropdims(maskNoneInf, dims=dim);
    maskInf = .!maskNoneInf;
    return maskInf
end

function getModelDistMatrix(data3d, latitudes)  
    # only take values where none (!) of the models has infinite values!! (Not just the two that are compared to one another)
    nbModels = size(data3d)[3];
    maskInf = getInfMask(data3d);

    # area weights differ across latitudes (because of the earth being a sphere)
    areaWeights = cos.(deg2rad.(latitudes));
    matrixS = zeros(nbModels, nbModels);
    # apply mask: ignore entries where at least one model has infinite value
    for i = 1:nbModels-1
        model_i = data3d[:,:,i];
        model_i[maskInf .== 1] .= 0;
        for j = (i + 1) : nbModels
            model_j = data3d[:, :, j];
            model_j[maskInf .== 1] .= 0;

            sqDiff = (model_i .- model_j).^2;
            values = areaWeights' .* sqDiff;
            s_ij = sqrt(sum(values));
            matrixS[i,j] = s_ij
        end
    end
    return matrixS
end


modelsTas, latitudes = loadModelData("tas");
modelsPr, latitudes = loadModelData("pr");
modelsPsl, latitudes = loadModelData("psl");

matrixTas = getModelDistMatrix(modelsTas, latitudes);
matrixPr = getModelDistMatrix(modelsPr, latitudes);
# make matrix symmetric
matrixTas = matrixTas .+ matrixTas';
matrixPr = matrixPr .+ matrixPr';

# Compare to processed data used in diagnostics
independenceTas = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_tas_CLIM.nc"), "dtas_CLIM");
independencePr = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_pr_CLIM.nc"), "dpr_CLIM");
precision = 4;
@assert round.(matrixTas, digits=precision) == round.(independenceTas, digits=precision)
@assert round.(matrixPr, digits=precision) == round.(independencePr, digits=precision)


# Data for overall mean plot
# These values are taken from the climwip recipe!
wPr, wTas = [0.25 0.5];
# normalize them 
wPr, wTas = [wPr wTas] ./ (wPr + wTas);

normalizedTas = matrixTas./median(matrixTas);
normalizedPr = matrixPr./median(matrixPr);
overallMeanIndependence = wPr .* normalizedPr .+ wTas .* normalizedTas;

independenceWeightsWork = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_overall_mean.nc"), "overall_mean");
@assert round.(overallMeanIndependence, digits=precision) == round.(independenceWeightsWork, digits=precision)


# 2.2. Performance weights 
# RMSE between models and observations

# get observational data from preproc-directory
function getObsData(climateVar)
    fn = "native6_ERA5_reanaly_v1_Amon_" * climateVar * "_1995-2014.nc";
    data = NetCDF.ncread(
        joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, climateVar * "_CLIM", fn),
        climateVar
    );
    latitudes = NetCDF.ncread(
        joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, climateVar * "_CLIM", fn),
        "lat"
    );
    return data, latitudes
end

# observational data
era5Tas, latitudes = getObsData("tas");
era5Psl, latitudes = getObsData("psl");
era5Pr, latitudes = getObsData("pr");


function getModelDataDist(models, observations, latitudes)      
    observations3d = reshape(observations, size(observations, 1), size(observations, 2), 1);
    areaWeights = cos.(deg2rad.(latitudes));

    # when using a single weightMatrix computed across all models and the observations
    #   weightMatrix = repeat(areaWeights', size(models, 1), 1);
    #   weightMatrix[maskInf .== 1] .= 0;

    distances = [];
    for i = 1 : size(models)[3]
        model_i = models[:,:,i];
        maskInf = getInfMask([observations3d;;;model_i]);
        maskedObs = deepcopy(observations)
        maskedObs[maskInf .== 1] .= 0;
    
        weightMatrix = repeat(areaWeights', size(models, 1), 1);
        weightMatrix[maskInf .== 1] .= 0;
        model_i[maskInf .== 1] .= 0;        
        sqDiff = (model_i .- maskedObs).^2;

        weightedValues = weightMatrix .* sqDiff;
        areaWeightedMean = sum(weightedValues) ./ sum(weightMatrix);
        distance = sqrt(areaWeightedMean);

        push!(distances, distance);
    end
    return distances
end

# Make sure to use a copy of the observational data, otherwise, it will be modified within the function by applying the mask!!
modelDataDistPr = getModelDataDist(modelsPr, deepcopy(era5Pr), latitudes);
modelDataDistTas = getModelDataDist(modelsTas, deepcopy(era5Tas), latitudes);
modelDataDistPsl = getModelDataDist(modelsPsl, deepcopy(era5Psl), latitudes);

# compare them to data from work dir
performancePr = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_pr_CLIM.nc"), "dpr_CLIM");
performanceTas = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_tas_CLIM.nc"), "dtas_CLIM");
performancePsl = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_psl_CLIM.nc"), "dpsl_CLIM");

@assert round.(performancePr, digits=precision) == round.(modelDataDistPr, digits=precision)
@assert round.(performanceTas, digits=precision) == round.(modelDataDistTas, digits=precision)
@assert round.(performancePsl, digits=precision) == round.(modelDataDistPsl, digits=precision)


# combine different performance diagnostics into an overall performance diagnostic:
# These values are taken from the climwip recipe!
wPr, wTas, wPsl = [2 1 1];

# normalize each diagnostic!
normalizedPr = performancePr./median(performancePr);
normalizedTas = performanceTas./median(performanceTas);
normalizedPsl = performancePsl ./ median(performancePsl);

# weight the normalized diagnostics according to weights from recipe and normalize
overallMeanPerformance = [wPr wTas wPsl] .* [normalizedPr normalizedTas normalizedPsl];
overallMeanPerformance = sum(overallMeanPerformance, dims=2) ./ sum([wPr wTas wPsl])

performanceWeightsWork = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_overall_mean.nc"), "overall_mean");
@assert vec(round.(overallMeanPerformance, digits=precision)) == round.(performanceWeightsWork, digits=precision)


# 2.3 Computation of final weights based on performance and independence weights
function averageEnsembleMembers(modelDiagnostics, modelNames)
    diagnosticDct = DefaultDict{String, Vector{Any}}(Vector{Any})
    for i in 1:length(modelNames)
        model = modelNames[i];
        diagnosticVal = modelDiagnostics[i];
        push!(diagnosticDct[model], diagnosticVal);
    end

    ensembleResults = [];
    uniqueModels = unique(modelNames);
    for model in uniqueModels
        push!(ensembleResults, mean(diagnosticDct[model]));
    end
    return ensembleResults, uniqueModels
end

modelNames = [split(model, "_")[2] for model in modelFileNames];
overallMeanPerformanceByEnsemble, models = averageEnsembleMembers(overallMeanPerformance, modelNames);


# TODO: this should be done with something similar to xarray with named dimensions!
# This is just a quiuck & dirty solution, knowing that the last 4 models are from the same ensemble
overallMeanIndependenceByEnsemble = overallMeanIndependence[1:3, 1:3];
ccsm4 = mean(overallMeanIndependence[1:3, 4:7], dims=2);
overallMeanIndependenceByEnsemble = hcat(overallMeanIndependenceByEnsemble, ccsm4);
overallMeanIndependenceByEnsemble = vcat(overallMeanIndependenceByEnsemble, [ccsm4' 0]);

# taken from climwip recipe!
sigmaD = 0.5;
sigmaS = 0.5;

numerator = exp.(-(overallMeanPerformanceByEnsemble./sigmaD).^2);
denominator = sum(exp.(-(overallMeanIndependenceByEnsemble./sigmaS).^2), dims=2); # (+1 in paper is from when model is compared to itself since exp(0)=1)

weights = numerator ./ denominator
weightsNormalized = weights ./ sum(weights)

# redistributed the weight for the 4 ensemble members
weightsComputed = vcat(weightsNormalized[1:3], fill(weightsNormalized[4]/4, 4));


weightsWork = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "weights.nc"), "weight")


@assert round.(weightsComputed, digits=precision) == round.(weightsWork, digits=precision)
