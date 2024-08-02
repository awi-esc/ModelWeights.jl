using SimilarityWeights
using NCDatasets
using DimensionalData
using Statistics

function checkPathToDir(path::String)
    return isdir(path)
end




function getData(pathToPreprocDir::String, climVar::String)
    if !checkPathToDir(pathToPreprocDir)
        throw(DomainError(pathToPreprocDir, "directory at given path does not exist!"))
    end 
    modelData =  SimilarityWeights.loadPreprocData(pathToPreprocDir, [climVar], "", ["CMIP6"]);
    data = modelData[climVar];
    return data
end

# precipitaton
#pathToPreprocDir = "/Users/brgrus001/output-from-albedo/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic"; 
pathToPreprocDir = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic";
data = getData(pathToPreprocDir, "pr");

means = dropdims(mean(data, dims=:model), dims=:model);
f1 = SimilarityWeights.plotMeansOnMap(means, "Precipitation means historical period")

# histogram of all data for specific location
# change longitudes to map from -180 to 180
longitudes = map(SimilarityWeights.lon360to180, Array(dims(data, :lon)));
data = set(data, :lon => longitudes);


# some locations
potsdam = Dict([("name", "Potsdam"), ("lat", 52.3906), ("lon", 13.0645)]);
veracruz = Dict([("name", "Veracruz"), ("lat", 19.1738), ("lon", -96.1342)]);
atacama = Dict([("name", "San Pedro de Atacama"), ("lat", -22.9087), ("lon", -68.1997)]);
mawsynram = Dict([("name", "Mawsynram"), ("lat", 25.300), ("lon", 91.583)]);

f2 = SimilarityWeights.plotHistAtPos(data, potsdam, "kg m-2 s-1")
f3 = SimilarityWeights.plotHistAtPos(data, veracruz, "kg m-2 s-1")
# in driest region of the earth
f4 = SimilarityWeights.plotHistAtPos(data, atacama, "kg m-2 s-1")
# appearently wettest place on earth
f5 = SimilarityWeights.plotHistAtPos(data, mawsynram, "kg m-2 s-1")


# AMOC (derived variable)
# pathToPreprocDir = "/Users/brgrus001/output-from-albedo/generated_recipe_historical_amoc_msftmz_20240730_090548/preproc/climatologic_diagnostic";
# pathToPreprocDir = "/Users/brgrus001/output-from-albedo/generated_recipe_historical_amoc_msftmz_20240731_153530/preproc/climatologic_diagnostic";
pathToPreprocDir = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/generated_recipe_historical_amoc_msftmz_20240731_153530/preproc/climatologic_diagnostic";

data = getData(pathToPreprocDir, "amoc");

SimilarityWeights.plotAMOC(data)


# AMOC (msftmz)


