using SimilarityWeights
using NCDatasets
using DimensionalData
using Statistics
using CairoMakie

function checkPathToDir(path::String)
    return isdir(path)
end


function getData(varToPath::Dict{String, String}, climVar::String, avgEnsembleMembers::Bool=true)
    modelData =  SimilarityWeights.loadPreprocData(varToPath, ["CMIP6"]);
    data = modelData[climVar];
    if avgEnsembleMembers
        data = SimilarityWeights.averageEnsembleVector(data, true)
    end
    return data
end

# precipitaton
varToPathPr =  Dict{String, String}("pr" => "/Users/brgrus001/output-from-albedo/historical/recipe_historical_pr/preproc/climatology_historical1/pr");
varToPathPr =  Dict{String, String}("pr" => "/Users/brgrus001/output-from-albedo/historical/recipe_historical_pr/preproc/climatology_historical3/pr");
output_dir = "/Users/brgrus001/output-from-albedo/" # TODO:adapt
#varToPath = Dict{String, String}("pr" => "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/generated_recipe_historical_pr_20240726_112338/preproc/climatologic_diagnostic");
data = getData(varToPathPr, "pr");

means = dropdims(mean(data, dims=:model), dims=:model);
f1 = SimilarityWeights.plotMeansOnMap(
    means, 
    "Precipitation means historical period", 
    SimilarityWeights.Target(
        directory  = output_dir;
        filename = "precipitation-historical1-unweighted-avg.png",
        save = true
    )
);

# histogram of all data for specific location
# change longitudes to map from -180 to 180
longitudes = SimilarityWeights.lon360to180.(Array(dims(data, :lon)));
data = set(data, :lon => longitudes);


# some locations
florida = Dict([("name", "Florida"), ("lat", 27.6648), ("lon", -81.5158)]);
potsdam = Dict([("name", "Potsdam"), ("lat", 52.3906), ("lon", 13.0645)]);
veracruz = Dict([("name", "Veracruz"), ("lat", 19.1738), ("lon", -96.1342)]);
# in driest region of the earth:
atacama = Dict([("name", "San Pedro de Atacama"), ("lat", -22.9087), ("lon", -68.1997)]);
# appearently wettest place on earth:
mawsynram = Dict([("name", "Mawsynram"), ("lat", 25.300), ("lon", 91.583)]);

f2 = SimilarityWeights.plotHistAtPos(data, potsdam)
f3 = SimilarityWeights.plotHistAtPos(data, veracruz)
f4 = SimilarityWeights.plotHistAtPos(data, atacama)
f5 = SimilarityWeights.plotHistAtPos(data, mawsynram)
f6 = SimilarityWeights.plotHistAtPos(data, florida)

# AMOC (derived variable)
# varToPath = Dict{String, String}("amoc" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_amoc_msftmz_20240730_090548/preproc/climatologic_diagnostic");
# varToPath = Dict{String, String}("amoc" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_amoc_msftmz_20240731_153530/preproc/climatologic_diagnostic");
# varToPath = Dict{String, String}("amoc" => "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/generated_recipe_historical_amoc_msftmz_20240731_153530/preproc/climatologic_diagnostic");

data = getData(varToPath, "amoc");
SimilarityWeights.convertKgsToSv!(data)
SimilarityWeights.plotAMOC(data)



# Models with different ensemble members

data_all_members = getData(varToPathPr, "pr", false);
longitudes = Array(dims(data_all_members, :lon));
if any(longitudes .> 180)
    longitudes = SimilarityWeights.lon360to180.(longitudes);
    data_all_members = set(data_all_members, :lon => SimilarityWeights.lon360to180.(Array(dims(data_all_members, :lon))));
end

loc = SimilarityWeights.getClosestGridPoint(veracruz, longitudes, Array(dims(data_all_members, :lat)));
loc = SimilarityWeights.getClosestGridPoint(potsdam, longitudes, Array(dims(data_all_members, :lat)));

SimilarityWeights.plotEnsembleSpread(data_all_members, loc["lon"], loc["lat"])