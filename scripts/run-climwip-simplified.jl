import ModelWeights as mw

using DimensionalData
using ColorSchemes
using NCDatasets
using Statistics

path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/climwip/climwip-simplified_20241013_073358"; 
path_configs = "/albedo/home/brgrus001/ModelWeights/configs/climwip_config";

# specify which data to load (default: all data will be loaded if compatible)
model_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var = false,
    is_model_data = true,
    subset = Dict("aliases" => ["calculate_weights_climwip"])
);

obs_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var = false,
    is_model_data = false,
    subset = Dict("aliases" => ["calculate_weights_climwip"])
);

# TODO: or load from seperate yaml file:
# data = mw.loadDataFromYAML(
# )

# Compute Weights
# configure parameters for computing weights
target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/";
fn_jld2 = "weights-climwip-simplified.jld2";
fn_nc = "weights-climwip-simplified.nc";
config_weights = mw.ConfigWeights(
    performance = Dict("tas_CLIM" => 1, "pr_CLIM" => 2, "psl_CLIM" => 1),
    independence = Dict("tas_CLIM" => 0.5, "pr_CLIM" => 0.25, "psl_CLIM" => 0),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    ref_perform_weights = "1995-2014",
    ref_indep_weights = "calculate_weights_climwip", # TODO: check this in paper!
    target_path = joinpath(target_dir, fn_jld2)
    # target_path = joinpath(target_dir, fn_nc)
);
weights = mw.computeWeights(
    collect(values(model_data)), collect(values(obs_data)), config_weights
);

# save weights as Julia object
mw.saveWeightsAsJuliaObj(weights, joinpath(target_dir, fn_jld2))
mw.saveWeightsAsNCFile(weights, joinpath(target_dir, fn_nc))

weights = mw.loadWeightsFromJLD2(joinpath(target_dir, fn_jld2));
weights_all_members = mw.distributeWeightsAcrossMembers(weights.w);

# make some Plots
figs_w = mw.plotWeights(weights.w; ylabel="overall weight")
figs_wP = mw.plotWeights(weights.wP; ylabel="Performance weight")
figs_wI = mw.plotWeights(weights.wI; ylabel="Independence weight")

mw.plotWeightContributions(weights.wI, weights.wP)

di_var = dropdims(
    reduce(+, weights.performance_distances, dims=:diagnostic), 
    dims=:diagnostic
);
figs_performance = mw.plotDistances(di_var, "Performance distances"; is_bar_plot = true);
figs_Di = mw.plotDistances(weights.Di, "Performance Di"; is_bar_plot = false);

figs_Sij = mw.plotDistancesIndependence(weights.Sij, dimname="model1");
sij_var = dropdims(reduce(+, weights.independence_distances, dims=:diagnostic), dims=:diagnostic)
figs_sij = mw.plotDistancesIndependence(sij_var, dimname="member1");


# Climwip Plots - Temperature map plots
# TODO fix:
# data_climwip = mw.loadDataFromESMValToolConfigs(
#     path_data, path_configs;
#     dir_per_var = false
# );

data_temp_map_future = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var=false,
    subset = Dict("aliases" => ["weighted_temperature_map_future"])
);
data_temp_map_reference = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var = false,
    subset = Dict("aliases" => ["weighted_temperature_map_reference"])
);
                
# compute weighted averages and plot results
data_ref = data_temp_map_reference["tas_CLIM_weighted_temperature_map_reference"].data;
data_future = data_temp_map_future["tas_CLIM_weighted_temperature_map_future"].data;
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
weighted_avg = mw.computeWeightedAvg(diff; weights = weights_all_members);
#weighted_avg = mw.sortLongitudesWest2East(weighted_avg);
f1 = mw.plotMeansOnMap(weighted_avg, title_f1; ColorSchemes.Reds.colors)

title_f2 = "Weighted minus unweighted mean temp. change: 2081-2100 minus 1995-2014";
unweighted_avg = mw.computeWeightedAvg(diff; use_members_equal_weights=false);
#unweighted_avg = mw.sortLongitudesWest2East(unweighted_avg);
diff_wu = weighted_avg .- unweighted_avg;
f2 = mw.plotMeansOnMap(diff_wu, title_f2; ColorSchemes.Reds.colors)


begin
    data_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_weighted_map.nc");
    data_ww_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(weighted_avg, data_ww_orig)
    # d = DimArray(data_ww_orig, (Dim{:lon}(Array(dims(weighted_avg, :lon))), Dim{:lat}(Array(dims(weighted_avg, :lat)))))
    # mw.plotMeansOnMap(d, title_f1 * " (original data)"; ColorSchemes.Reds.colors);
end

begin
    data_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_difference_map.nc");
    data_wu_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(diff_wu, data_wu_orig)
end


# Apply computed weights - Temperature graph plots
data_temp_graph = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var = false,
    subset = Dict("aliases" => ["weighted_temperature_graph"])
);
data_graph = data_temp_graph["tas_ANOM_weighted_temperature_graph"].data;
# this will compute the weighted avg based on the average across the respective members of each model
#weighted_avg = mw.applyWeights(data_graph, weights_all_members);

weighted_avg = mw.computeWeightedAvg(data_graph; weights = weights_all_members);
unweighted_avg = mw.computeWeightedAvg(data_graph; use_members_equal_weights = false);
uncertainties = mw.getUncertaintyRanges(data_graph, weights_all_members);

f3 = mw.plotTempGraph(
    data_graph, 
    (weighted=weighted_avg, unweighted=unweighted_avg),
    uncertainties,
    "Temperature anomaly relative to 1981-2010";
    ylabel = "Temperature anomaly"
)
tas_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_graph/weighted_temperature_graph/temperature_anomalies.nc")["tas"];
@assert tas_orig == data_graph

unc_unweighted_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/orig-data-temp-graph/uncertainty_range.nc");
compareToOrigData(unc_unweighted_orig["tas"][:,:][1,:], map(x -> x[1], uncertainties.unweighted))
compareToOrigData(unc_unweighted_orig["tas"][:,:][2,:], map(x -> x[2], uncertainties.unweighted))
unc_weighted_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/orig-data-temp-graph/uncertainty_range_weighted.nc");
compareToOrigData(unc_weighted_orig["tas"][:,:][1,:], map(x -> x[1], uncertainties.weighted))
# TODO: only this is not equal for a handful of indices!
uncertainties_weighted_orig = unc_weighted_orig["tas"][:,:][2,:];
uncertainties_weighted = map(x -> x[2], uncertainties.weighted);
diff = uncertainties_weighted .- uncertainties_weighted_orig;
indices = findall(x -> x>0.0001, diff)
diff[indices]
compareToOrigData(uncertainties_weighted_orig, uncertainties_weighted)

# data_graph_models = mw.summarizeEnsembleMembersVector(data_graph, true);
# uncertainties_models = mw.getUncertaintyRanges(data_graph_models, weights.w);


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
