using NetCDF
using Statistics

include("load-data-utils.jl");

# 1. Preprocessed data for weighted_temperature_graph diagnostics
path2preprocWeightedTempGraph = joinpath(PATH_TO_PREPROC_DIR, "weighted_temperature_graph", "tas")
preprocTempGraph, _ = loadNCdataInDir(path2preprocWeightedTempGraph, "tas", []);

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
# Model-Model Distance matrices 
PATH_TO_PREPROC_WEIGHTS_DIR = joinpath(PATH_TO_PREPROC_DIR, "calculate_weights_climwip");

preprocTas, filenamesTas = loadNCdataInDir(joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, "tas_CLIM"), "tas", ["CMIP"], false);
preprocPr, filenamesPr = loadNCdataInDir(joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, "pr_CLIM"), "pr", ["CMIP"], false);

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



function getModelDistMatrix(climateVar::String)  
    climateVarDir = climateVar * "_CLIM";
    preprocData, filenames = loadNCdataInDir(joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, climateVarDir), climateVar, ["CMIP"], false);
    latitudes = NetCDF.ncread(joinpath(PATH_TO_PREPROC_WEIGHTS_DIR, climateVarDir, filenames[1]), "lat");

    # only take values where none (!) of the models has infinite values!! (Not just the two that are compared to one another)
    data = cat(preprocData..., dims=3);
    nbModels = size(data)[3];
    maskNotInf = data .!= 1.0f20;
    maskNoneInf = sum(maskNotInf, dims=3) .== nbModels;  
    maskNoneInf = dropdims(maskNoneInf, dims=3);
    maskInf = .!maskNoneInf;

    # area weights differ across latitudes (because of the earth being a sphere)
    areaWeights = cos.(deg2rad.(latitudes));
    matrixS = zeros(nbModels, nbModels);
    # apply mask: ignore entries where one model has infinite value
    for i = 1:nbModels-1
        model_i = data[:,:,i];
        model_i[maskInf .== 1] .= 0;
        for j = (i + 1) : nbModels
            model_j = data[:, :, j];
            model_j[maskInf .== 1] .= 0;

            sqDiff = (model_i .- model_j).^2;
            values = areaWeights' .* sqDiff;
            s_ij = sqrt(sum(values));
            matrixS[i,j] = s_ij
        end
    end
    return matrixS
end


# make matrix symmetric
matrixTas = getModelDistMatrix("tas");
matrixPr = getModelDistMatrix("pr");

matrixTas = matrixTas .+ matrixTas';
matrixPr = matrixPr .+ matrixPr';

# Compare to processed data used in diagnostics
independenceTas = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_tas_CLIM.nc"), "dtas_CLIM");
independencePr = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_pr_CLIM.nc"), "dpr_CLIM");


precision = 4;
@assert round.(matrixTas, digits=precision) == round.(independenceTas, digits=precision)
@assert round.(matrixPr, digits=precision) == round.(independencePr, digits=precision)







