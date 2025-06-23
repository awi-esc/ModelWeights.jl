import ModelWeights.Data as mwd
import ModelWeights.Weights as mww
import ModelWeights.Plots as mwp

using CairoMakie
using Colors
using ColorSchemes
using DimensionalData
using NCDatasets

# Load data for computing weights (output from ESMValTool recipe)
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/climwip/climwip-simplified_20241013_073358"; 
path_configs = "./configs/climwip_config";

data = mwd.loadDataFromESMValToolRecipes(
    path_data, path_configs;
    dir_per_var = false,
    subset = Dict("aliases" => ["calculate_weights_climwip"])
)

# Compute Weights (performance based on historical period)
target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/";
fn_jld2 = "weights-climwip-simplified.jld2";
fn_nc = "weights-climwip-simplified.nc";
config_weights = mww.ConfigWeights(
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

dists_perform = mwd.distancesData(data.models, data.obs, config_weights.performance);
dists_indep = mwd.distancesModels(data.models, config_weights.independence);
weights = mww.climwipWeights(dists_indep, dists_perform, config_weights);


# Load weights from disk
weights = mw.readDataFromDisk(config_weights.target_path; variable="weights");

# make some Plots
fig_weights = mwp.plotWeights(weights.w; title="Climwip test basic; weights")
# mwp.plotWeightContributions(weights.wI, weights.wP)
# figs_performance = mwp.plotDistancesByVar(
#     weights.performance_distances, "Performance distances "; is_bar_plot=true
# )

figs_Sij = mwp.plotDistancesIndependence(weights.Sij, "model1")
figs_Sij[1]
sij_var = dropdims(
    reduce(+, weights.independence_distances, dims=:diagnostic), 
    dims=:diagnostic
)
# figs_sij = mw.plotDistancesIndependence(sij_var, "member1")


# Climwip Plots - Temperature map plots
data_temp_map_future = mwd.loadDataFromESMValToolRecipes(
    path_data, path_configs;
    dir_per_var=false,
    subset = Dict("aliases" => ["weighted_temperature_map_future"])
)
data_temp_map_reference = mwd.loadDataFromESMValToolRecipes(
    path_data, path_configs;
    dir_per_var = false,
    subset = Dict("aliases" => ["weighted_temperature_map_reference"])
)
                
# compute weighted averages and plot results
data_ref = data_temp_map_reference.models["tas_CLIM_weighted_temperature_map_reference"];
data_future = data_temp_map_future.models["tas_CLIM_weighted_temperature_map_future"];
# just to align with original data
data_ref = data_ref[lat = Where(x -> x <= 68.75)];
data_future = data_future[lat = Where(x -> x <= 68.75)];
diff = data_future .- data_ref;

# sanity checks: compare to original data (must have adapted the latitudes above)
function compareToOrigData(data, data_orig)
    mat = isapprox.(data, data_orig, atol=10^-4);
    @assert ismissing.(data_orig) == ismissing.(data)
    approxEqual = all(x -> ismissing(x) || x, mat)
    @assert approxEqual
end

title_f1 = "Weighted mean temp. change 2081-2100 minus 1995-2014";
weighted_avg = mww.weightedAvg(diff; weights = weights.w_members);
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
f1

title_f2 = "Weighted minus unweighted mean temp. change: 2081-2100 minus 1995-2014";
unweighted_avg = mww.weightedAvg(diff; use_members_equal_weights=false);
#unweighted_avg = mw.sortLongitudesWest2East(unweighted_avg);
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
f2

begin
    data_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_weighted_map.nc");
    data_ww_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(weighted_avg, data_ww_orig)
    # d = DimArray(data_ww_orig, (Dim{:lon}(Array(dims(weighted_avg, :lon))), Dim{:lat}(Array(dims(weighted_avg, :lat)))))
    # f = Figure();
    # mw.plotValsOnMap!(f, d, title_f1 * " (original data)"; ColorSchemes.Reds.colors);
end

begin
    data_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_difference_map.nc");
    data_wu_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(diff_wu, data_wu_orig)
end


# Apply computed weights - Temperature graph plots
data_temp_graph = mwd.loadDataFromESMValToolRecipes(
    path_data, path_configs;
    dir_per_var = false,
    subset = Dict("aliases" => ["weighted_temperature_graph"])
);
data_graph = data_temp_graph.models["tas_ANOM_weighted_temperature_graph"];
# this will compute the weighted avg based on the average across the respective members of each model
#weighted_avg = mw.applyWeights(data_graph, weights.w_members);

weighted_avg = mww.weightedAvg(data_graph; weights = weights.w_members);
unweighted_avg = mww.weightedAvg(data_graph; use_members_equal_weights = false);

uncertainties_weighted = mwd.uncertaintyRanges(data_graph; w = weights.w_members);
uncertainties_unweighted = mwd.uncertaintyRanges(data_graph);

f3 = mw.plotTempGraph(
    data_graph, 
    (weighted=weighted_avg, unweighted=unweighted_avg),
    (weighted=uncertainties_weighted, unweighted=uncertainties_unweighted),
    "Temperature anomaly relative to 1981-2010";
    ylabel = "Temperature anomaly"
)
tas_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_graph/weighted_temperature_graph/temperature_anomalies.nc")["tas"];
@assert tas_orig == data_graph

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

data_graph_models = mwd.summarizeEnsembleMembersVector(data_graph, true);

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
