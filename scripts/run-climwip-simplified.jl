import ModelWeights as mw
import ModelWeights.Data as mwd
import ModelWeights.Weights as mww
import ModelWeights.Plots as mwp

using CairoMakie
using Colors
using ColorSchemes
using DimensionalData
using NCDatasets
using NetCDF
using Statistics
using YAXArrays

# Load data for computing weights (output from ESMValTool recipe)
path_data = "test/data/climwip-simplified"; 
path_recipes = "./configs/climwip_config";

plot_dir = "reproduce-climwip-figs"

model_data = mw.defineDataMap(
    path_data, 
    path_recipes, 
    :esmvaltool_recipes; 
    dir_per_var = false,
    dtype = "cmip",
    filename_format = :esmvaltool,
    constraint = Dict(
        "aliases" => ["calculate_weights_climwip"], 
        "mips" => ["CMIP5", "CMIP6"]
    )
)
obs_data =  mw.defineDataMap(
    path_data, 
    path_recipes, 
    :esmvaltool_recipes; 
    dir_per_var = false,
    filename_format = :esmvaltool,
    dtype = "obs",
    constraint = Dict(
        "aliases" => ["calculate_weights_climwip"],
        "filenames" => ["ERA5"]
    )
)
mwd.apply!(obs_data, mwd.setDim, :model, nothing, ["ERA5"])


# Compute Weights (performance based on historical period)
config = mww.ConfigWeights(
    performance = Dict(
        "tas_CLIM_calculate_weights_climwip" => 1, 
        "pr_CLIM_calculate_weights_climwip" => 2, 
        "psl_CLIM_calculate_weights_climwip" => 1
    ),
    independence = Dict(
        "tas_CLIM_calculate_weights_climwip" => 0.5,
        "pr_CLIM_calculate_weights_climwip" => 0.25,
        "psl_CLIM_calculate_weights_climwip" => 0
    ),
    sigma_independence = 0.5,
    sigma_performance = 0.5
);
# this is done just to compare with their data which is sorted that way
mwd.apply!(model_data, mwd.sortLongitudesEast2West)
mwd.apply!(obs_data, mwd.sortLongitudesEast2West)
weights = mww.climwipWeights(model_data, obs_data, config; suffix_name="simplified", norm_avg_members = false);


# check performance distances
PATH_TO_WORK_DIR = "reproduce-climwip-figs/recipe_climwip_test_basic_data/work"
models = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_tas_CLIM.nc"), "model_ensemble")
performanceTas = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_tas_CLIM.nc"), "dtas_CLIM");
performancePr = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_pr_CLIM.nc"), "dpr_CLIM");
performancePsl = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_psl_CLIM.nc"), "dpsl_CLIM");

dTas = weights.performance_distances["tas_CLIM_calculate_weights_climwip"].data
dPr = weights.performance_distances["pr_CLIM_calculate_weights_climwip"].data
dPsl = weights.performance_distances["psl_CLIM_calculate_weights_climwip"].data
all(isapprox.(performanceTas, dTas))
all(isapprox.(performancePr, dPr))
all(isapprox.(performancePsl, dPsl))

dPerformance = [dTas, dPr, dPsl] ./ median.([dTas, dPr, dPsl])
Di = sum([1/4, 2/4, 1/4] .* dPerformance)
performanceOverall = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip/climwip/performance_overall_mean.nc"), "overall_mean")

all(isapprox.(Di, performanceOverall))
# our computed generalized distances yield slightly different results (if norm_avg_members was true):
weights.Di.data 
# this is because we first summarize members into one average distance for each model and then take the median
# instead of taking the median over all models and members:

# check independence distances
independenceTas = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip/climwip/independence_tas_CLIM.nc"), "dtas_CLIM")
independencePr = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip/climwip/independence_pr_CLIM.nc"), "dpr_CLIM")

dTas = weights.independence_distances["tas_CLIM_calculate_weights_climwip"].data
dPr = weights.independence_distances["pr_CLIM_calculate_weights_climwip"].data
# the following is different because we normalize the area weights, which they don't for the independence weights (using pdist with unnormalized weights)
# (if we don't normalize the weights we get the same results as they do)
all(isapprox.(independenceTas, dTas))
all(isapprox.(independencePr, dPr))

# yet the computed generalized distances are identical because the normalization of the median cancels out, since our distance values and their distance values
# just differ by a factor (the same across all model pairs):
independenceTas ./ dTas
independencePr ./ dPr

dIndependence = [dTas, dPr] ./ median.([dTas, dPr])
Sij = sum([2/3, 1/3] .* dIndependence)
independence_overall = NetCDF.ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip/climwip/independence_overall_mean.nc"), "overall_mean")
all(isapprox.(Sij, independence_overall))
# our generalized distances are again slightly different because we first summarize members into one average distance for each model and then take the median
# instead of taking the median over all models and members (as we did here):
weights.Sij.data


# results from ESMValTool run
orig_wP = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip/climwip/performance_overall_mean.nc"))
orig_wI = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip/climwip/independence_overall_mean.nc"))
orig_wIP = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip/climwip/weights.nc"))
models = orig_wIP["model_ensemble"][:]
orig_w = YAXArray((Dim{:member}(models),), orig_wIP["weight"][:])
orig_w.data
# comparison with our weights (slightly different for reasons mentioned above)
weights.w[weight = At("wIP-simplified")].data

# make some Plots
fig_weights = mwp.plotWeights(weights.w; title="Climwip test basic; weights")
mwp.savePlot(fig_weights, joinpath(plot_dir, "weights.png"); overwrite=true)
# plot generalized distances for independence
figs_Sij = mwp.plotDistancesIndependence(weights.Sij, "model1")
mwp.savePlot(figs_Sij[1], joinpath(plot_dir, "Sij.png"); overwrite=true)
# plot Moel-Model distances for all members for each diagnostic
dists = weights.independence_distances["tas_CLIM_calculate_weights_climwip"]
figs_sij = mwp.plotDistancesIndependence(dists, "member1"; title="Model-Model distances for tas")
mwp.savePlot(figs_sij[1], joinpath(plot_dir, "sij-members-tas.png"); overwrite=true)

dists = weights.independence_distances["pr_CLIM_calculate_weights_climwip"]
figs_sij = mwp.plotDistancesIndependence(dists, "member1"; title="Model-Model distances for pr")
mwp.savePlot(figs_sij[1], joinpath(plot_dir, "sij-members-pr.png"); overwrite=true)

# plot generalized distances for performance
fig_Di = mwp.plotDistances(weights.Di, "Generalized distances Di")
mwp.savePlot(fig_Di, joinpath(plot_dir, "Di.png"); overwrite=true)
# plot distances for all members for each diagnostic
fig_di = mwp.plotDistances(weights.performance_distances["tas_CLIM_calculate_weights_climwip"], "Model-data distances for tas")
mwp.savePlot(fig_di, joinpath(plot_dir, "di-members-tas.png"); overwrite=true)

fig_di = mwp.plotDistances(weights.performance_distances["pr_CLIM_calculate_weights_climwip"], "Model-data distances for pr")
mwp.savePlot(fig_di, joinpath(plot_dir, "di-members-pr.png"); overwrite=true)

fig_di = mwp.plotDistances(weights.performance_distances["psl_CLIM_calculate_weights_climwip"], "Model-data distances for psl")
mwp.savePlot(fig_di, joinpath(plot_dir, "di-members-psl.png"); overwrite=true)


# Climwip Plots - Temperature map plots
data_temp_map_future = mw.defineDataMap(
    path_data, path_recipes, :esmvaltool_recipes; 
    dtype = "cmip",
    dir_per_var = false,
    constraint = Dict(
        "aliases" => ["weighted_temperature_map_future"],
        "mips" => ["CMIP5", "CMIP6"]
    )
)
data_temp_map_reference = mw.defineDataMap(
    path_data, path_recipes, :esmvaltool_recipes;
    dir_per_var = false,
    dtype = "cmip",
    constraint = Dict("aliases" => ["weighted_temperature_map_reference"])
)
              
# compute weighted averages and plot results
data_ref = data_temp_map_reference["tas_CLIM_weighted_temperature_map_reference"];
data_future = data_temp_map_future["tas_CLIM_weighted_temperature_map_future"];
# just to align with original data
data_ref = data_ref[lat = Where(x -> x <= 68.75)];
data_future = data_future[lat = Where(x -> x <= 68.75)];

diff = data_future.data .- data_ref.data;
members = mwd.sharedModels(Dict{String, YAXArray}(weights.performance_distances), :member)
w_members = mww.distributeWeightsAcrossMembers(weights.w[weight = At("wIP-simplified")], members)

title_f1 = "Weighted mean temp. change 2081-2100 minus 1995-2014";
weighted_avg = mww.weightedAvg(YAXArray(dims(diff), collect(diff)); weights = w_members);
#weighted_avg = mw.sortLongitudesWest2East(weighted_avg);
f1 = Figure();
cmap = cgrad([:white, :salmon, :red, :darkred], 10; categorical = true)
mwp.plotValsOnMap!(
    f1, weighted_avg, title_f1; 
    color_range = (2.5, 6.5),
    colors = cmap[2:end-1],
    high_clip = cmap[end],
    low_clip = cmap[1],
    xlabel_rotate = false,
    east_west_labels = true
)
mwp.savePlot(f1, joinpath(plot_dir, "weighted-mean-temp-change.png"))

title_f2 = "Weighted minus unweighted mean temp. change: 2081-2100 minus 1995-2014";
unweighted_avg = mww.weightedAvg(diff; use_members_equal_weights=false);
#unweighted_avg = mwd.sortLongitudesWest2East(unweighted_avg);
diff_wu = weighted_avg .- unweighted_avg;
f2 = Figure();
cmap = reverse(Colors.colormap("RdBu", logscale=false, mid=0.5))
mwp.plotValsOnMap!(f2, diff_wu, title_f2; 
    color_range = (-1, 1),
    colors = cmap[2:end-1],
    high_clip = cmap[end],
    low_clip = cmap[1],
    xlabel_rotate = false,
    east_west_labels = true
)
mwp.savePlot(f2, joinpath(plot_dir, "weighted-minus-unweighted-temp-change.png"))

# sanity checks: compare to original data (must have adapted the latitudes above)
function compareToOrigData(data, data_orig)
    mat = isapprox.(data, data_orig, atol=10^-4);
    @assert ismissing.(data_orig) == ismissing.(data)
    approxEqual = all(x -> ismissing(x) || x, mat)
    @assert approxEqual
end

begin
    data_orig = NCDataset(joinpath(PATH_TO_WORK_DIR, "weighted_temperature_map/weighted_temperature_map/temperature_change_weighted_map.nc"));
    data_ww_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(Array(weighted_avg), data_ww_orig)
end

begin
    data_orig = NCDataset(joinpath(PATH_TO_WORK_DIR, "weighted_temperature_map/weighted_temperature_map/temperature_change_difference_map.nc"));
    data_wu_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(Array(diff_wu), data_wu_orig)
end


# Apply computed weights - Temperature graph plots
data_temp_graph = mw.defineDataMap(
    path_data, 
    path_recipes,
    :esmvaltool_recipes;
    dir_per_var = false,
    constraint = Dict("aliases" => ["weighted_temperature_graph"])
);
# transpose brings time to second dimension, otherwise slicing in uncertaintyRanges fails..
data_graph = transpose(data_temp_graph["tas_ANOM_weighted_temperature_graph"]);
# this will compute the weighted avg based on the average across the respective members of each model
#weighted_avg = mww.applyWeights(data_graph, w_members);

weighted_avg = mww.weightedAvg(data_graph; weights = w_members);
unweighted_avg = mww.weightedAvg(data_graph; use_members_equal_weights = false);

uncertainties_weighted = mwd.uncertaintyRanges(data_graph; w = w_members);
uncertainties_unweighted = mwd.uncertaintyRanges(data_graph);

f3 = mwp.plotTempGraph(
    data_graph, 
    (weighted=weighted_avg, unweighted=unweighted_avg),
    (weighted=uncertainties_weighted, unweighted=uncertainties_unweighted),
    "Temperature anomaly relative to 1981-2010";
    ylabel = "Temperature anomaly"
)
mwp.savePlot(f3, joinpath(plot_dir, "temp-graph-anomalies.png"))
tas_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_graph/weighted_temperature_graph/temperature_anomalies.nc")["tas"];
@assert Array(tas_orig["tas"]) == Array(transpose(data_graph))

# Recheck the following:
#unc_unweighted_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/orig-data-temp-graph/uncertainty_range.nc");
# compareToOrigData(unc_unweighted_orig["tas"][:,:][1,:], map(x -> x[1], uncertainties_unweighted))
# compareToOrigData(unc_unweighted_orig["tas"][:,:][2,:], map(x -> x[2], uncertainties_unweighted))
# unc_weighted_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/orig-data-temp-graph/uncertainty_range_weighted.nc");
# compareToOrigData(unc_weighted_orig["tas"][:,:][1,:], map(x -> x[1], uncertainties_weighted))
# # # TODO: only this is not equal for a handful of indices!
# uncertainties_weighted_orig = unc_weighted_orig["tas"][:,:][2,:];
# uncertainties_weighted = map(x -> x[2], uncertainties_weighted);
# diff = uncertainties_weighted .- uncertainties_weighted_orig;
# indices = findall(x -> x>0.0001, diff)
# diff[indices]
# compareToOrigData(uncertainties_weighted_orig, uncertainties_weighted)

data_graph_models = mwd.summarizeMembers(data_graph);

weighted_avg_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/orig-data-temp-graph/central_estimate_weighted.nc")
compareToOrigData(weighted_avg_orig["tas"][:], weighted_avg[:])
unweighted_avg_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/orig-data-temp-graph/central_estimate.nc")
compareToOrigData(unweighted_avg_orig["tas"][:], unweighted_avg[:])




# TODO: Plot ensemble spread for some models with >1 ensemble member
# data = modelDataRef["CLIM"]["tas"]
# models_kept = ["EC-Earth3", "ACCESS-CM2"]
# indices = findall(m -> m in models_kept, Array(dims(data, :model)))
# shared_models = data.metadata["full_model_names"][indices]
# data = sw.getModelSubset(data, shared_models)
# sw.plotEnsembleSpread(data, 7.5, 82.5)
