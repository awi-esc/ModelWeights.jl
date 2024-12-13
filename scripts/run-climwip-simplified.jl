import ModelWeights as mw

using DimensionalData
using ColorSchemes
using NCDatasets

path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/climwip/climwip-simplified_20241013_073358"; 
path_configs = "/albedo/home/brgrus001/ModelWeights/configs/climwip_config";

# specify which data to load (default: all data will be loaded if compatible)
model_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var = false,
    is_model_data = true,
    only_shared_models = false,
    subset = mw.Constraint(aliases = ["calculate_weights_climwip"])
);

obs_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var = false,
    is_model_data = false,
    subset = mw.Constraint(aliases = ["calculate_weights_climwip"])
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
    ref_period = "1995-2014", 
    target_path = joinpath(target_dir, fn_jld2)
    # target_path = joinpath(target_dir, fn_nc)
);
weights = mw.computeWeights(model_data, obs_data, config_weights);

# save weights as Julia object
mw.saveWeightsAsJuliaObj(weights, joinpath(target_dir, fn_jld2))
mw.saveWeightsAsNCFile(weights, joinpath(target_dir, fn_nc))
weights = mw.loadWeightsFromJLD2(target_path)
weights_all_members = mw.distributeWeightsAcrossMembers(weights.w);

# make some Plots
figs_w = mw.plotWeights(weights.w; label="overall weight")
figs_wP = mw.plotWeights(weights.wP; label="Performance weight")
figs_wI = mw.plotWeights(weights.wI; label="Independence weight")

mw.plotWeightContributions(weights.wI, weights.wP)

di_var = dropdims(
    reduce(+, weights.performance_distances, dims=:diagnostic), 
    dims=:diagnostic
);
figs_performance = mw.plotDistancesPerformance(di_var; is_bar_plot = true);
figs_Di = mw.plotDistancesPerformance(weights.Di; is_bar_plot = false);

figs_Sij = mw.plotDistancesIndependence(weights.Sij, dimname="model1");
sij_var = dropdims(reduce(+, weights.independence_distances, dims=:diagnostic), dims=:diagnostic)
figs_sij = mw.plotDistancesIndependence(sij_var, dimname="member1");


# Apply computed weights - Temperature map plots
# TODO fix:
# data_climwip = mw.loadDataFromESMValToolConfigs(
#     path_data, path_configs;
#     dir_per_var = false
# );

data_temp_map_future = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var=false,
    subset = mw.Constraint(aliases = ["weighted_temperature_map_future"])
);
data_temp_map_reference = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var = false,
    only_shared_models = true,
    subset = mw.Constraint(aliases = ["weighted_temperature_map_reference"])
);
                
# compute weighted averages and plot results
data_ref = mw.indexData(data_temp_map_reference, "tas", "CLIM", "weighted_temperature_map_reference")
data_future = mw.indexData(data_temp_map_future, "tas", "CLIM", "weighted_temperature_map_future")
# to align with original data
data_ref = data_ref[lat = Where(x -> x <= 68.75)];
data_future = data_future[lat = Where(x -> x <= 68.75)];

weighted_ref = mw.computeWeightedAvg(data_ref; weights = weights_all_members)
weighted_future = mw.computeWeightedAvg(data_future; weights = weights_all_members)
diff_ww = weighted_future .- weighted_ref;
diff_ww = mw.sortLongitudesWest2East(diff_ww);

# Weighted mean temp. change 2081-2100 minus 1995-2014
f1 = mw.plotMeansOnMap(
    diff_ww, "Weighted mean temp. change 2081-2100 minus 1995-2014";
    ColorSchemes.Reds.colors
)

# some sanity checks
function compareToOrigData(data, data_orig)
    mat = isapprox.(data, data_orig, atol=10^-2);
    @assert ismissing.(data_orig) == ismissing.(data)
    approxEqual = all(x -> ismissing(x) || x, mat)
    @assert approxEqual
end

begin
    # compare to original data (must have adapted the latitudes above)
    data_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_weighted_map.nc");
    data_ww_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(diff_ww, data_ww_orig)
    d = DimArray(data_ww_orig, (Dim{:lon}(Array(dims(diff_ww, :lon))), Dim{:lat}(Array(dims(diff_ww, :lat)))))
    mw.plotMeansOnMap(d, "Weighted mean temp. change 2081-2100 - 1995-2014"; ColorSchemes.Reds.colors);
end

####### DOESNT WORK YET ########
# Weighted minus unweighted mean temp. change 2081-2100 minus 1995-2014
data_ref_models = mw.summarizeEnsembleMembersVector(data_ref, true);
data_future_models = mw.summarizeEnsembleMembersVector(data_future, true);

unweighted_ref = mw.computeWeightedAvg(data_ref);
unweighted_future = mw.computeWeightedAvg(data_future);
diff_uu = unweighted_future .- unweighted_ref;
diff_uu = mw.sortLongitudesWest2East(diff_uu);
diff_wu = diff_ww .- diff_uu;

f2 = mw.plotMeansOnMap(
    diff_wu,
    "Weighted minus unweighted mean temp. change: 2081-2100 minus 1995-2014";
    ColorSchemes.Reds.colors
)

# compare to orig data
begin
    data_orig = NCDataset("/albedo/home/brgrus001/ModelWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_difference_map.nc");
    data_wu_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(diff_wu, data_wu_orig)
end
####### DOESNT WORK YET: results not equal (AssertionError) ########


# Apply computed weights - Temperature graph plots
data_temp_graph = mw.loadDataFromESMValToolConfigs(
    path_data, path_configs;
    dir_per_var = false,
    subset = mw.Constraint(aliases = ["weighted_temperature_graph"])
);
data_graph = mw.indexData(data_temp_graph, "tas", "ANOM", "weighted_temperature_graph")
# TODO: check applyWeights fn and difference!
# this will compute the weighted avg based on the average across the respective members of each model
#weighted_avg = mw.applyWeights(data_graph, weights_all_members);

weighted_avg = mw.computeWeightedAvg(data_graph; weights = weights_all_members);
unweighted_avg = mw.computeWeightedAvg(data_graph);

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

# TODO check if uncertainties are also identical! And weighted/unweighted avgs!


# TODO: Plot ensemble spread for some models with >1 ensemble member
# data = modelDataRef["CLIM"]["tas"]
# models_kept = ["EC-Earth3", "ACCESS-CM2"]
# indices = findall(m -> m in models_kept, Array(dims(data, :model)))
# shared_models = data.metadata["full_model_names"][indices]
# data = sw.getModelSubset(data, shared_models)
# sw.plotEnsembleSpread(data, 7.5, 82.5)
