import SimilarityWeights as sw
using NCDatasets
using DimensionalData
using Statistics
using CairoMakie

function checkPathToDir(path::String)
    return isdir(path)
end


function getData(
    varToPath::Dict{String, Dict{String, String}}, 
    diagnostic::String, 
    climVar::String, 
    avgEnsembleMembers::Bool=true
)
    modelData =  sw.loadPreprocData(varToPath; included=["CMIP6"]);
    data = modelData[diagnostic][climVar];
    if avgEnsembleMembers
        data = sw.averageEnsembleVector(data, true)
    end
    return data
end

# precipitaton
basePath = "/Users/brgrus001/output-from-albedo/historical/recipe_historical_pr"
period = "historical1"
climVar = "pr"
diagnostic = "CLIM"
varToPath = Dict(diagnostic => Dict(climVar => joinpath(basePath, "preproc", period, climVar * "_" * diagnostic)))
data = getData(varToPath, diagnostic, climVar);

means = dropdims(mean(data, dims=:model), dims=:model);
f1 = sw.plotMeansOnMap(
    means, 
    "Precipitation means historical period"
);
sw.savePlot(f1, "output", "precipitation-historical1-unweighted-avg.png")


# histogram of all data for specific location
# change longitudes to map from -180 to 180
longitudes = sw.lon360to180.(Array(dims(data, :lon)));
data = set(data, :lon => longitudes);


# some locations
florida = Dict([("name", "Florida"), ("lat", 27.6648), ("lon", -81.5158)]);
potsdam = Dict([("name", "Potsdam"), ("lat", 52.3906), ("lon", 13.0645)]);
veracruz = Dict([("name", "Veracruz"), ("lat", 19.1738), ("lon", -96.1342)]);
# in driest region of the earth:
atacama = Dict([("name", "San Pedro de Atacama"), ("lat", -22.9087), ("lon", -68.1997)]);
# appearently wettest place on earth:
mawsynram = Dict([("name", "Mawsynram"), ("lat", 25.300), ("lon", 91.583)]);

f2 = sw.plotHistAtPos(data, potsdam)
f3 = sw.plotHistAtPos(data, veracruz)
f4 = sw.plotHistAtPos(data, atacama)
f5 = sw.plotHistAtPos(data, mawsynram)
f6 = sw.plotHistAtPos(data, florida)


###############################################################################
# TODO: update AMOC Plot
# AMOC (derived variable)
# varToPath = Dict{String, String}("amoc" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_amoc_msftmz_20240730_090548/preproc/climatologic_diagnostic");
# varToPath = Dict{String, String}("amoc" => "/Users/brgrus001/output-from-albedo/generated_recipe_historical_amoc_msftmz_20240731_153530/preproc/climatologic_diagnostic");
# varToPath = Dict{String, String}("amoc" => "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/generated_recipe_historical_amoc_msftmz_20240731_153530/preproc/climatologic_diagnostic");

# data = getData(varToPath, "amoc");
# sw.convertKgsToSv!(data)
# sw.plotAMOC(data)
###############################################################################


# Models with different ensemble members
data_all_members = getData(varToPath, "CLIM", "pr", false);
longitudes = Array(dims(data_all_members, :lon));
if any(longitudes .> 180)
    longitudes = sw.lon360to180.(longitudes);
    data_all_members = set(data_all_members, :lon => sw.lon360to180.(Array(dims(data_all_members, :lon))));
end

loc = sw.getClosestGridPoint(veracruz, longitudes, Array(dims(data_all_members, :lat)));
loc = sw.getClosestGridPoint(potsdam, longitudes, Array(dims(data_all_members, :lat)));

sw.plotEnsembleSpread(data_all_members, loc["lon"], loc["lat"])