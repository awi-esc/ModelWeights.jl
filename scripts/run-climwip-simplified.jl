import SimilarityWeights as sw

using DimensionalData
using ColorSchemes
using NCDatasets

base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/climwip/climwip-simplified_20241013_073358"; 
config_path = "/albedo/home/brgrus001/SimilarityWeights/configs/climwip_config";

# specify which data to load (default: all data will be loaded if compatible)
data_weights = sw.loadData(
    base_path,
    config_path;
    dir_per_var=false,
    common_models_across_vars=true,
    subset=Dict(
        #"variables" => Vector{String}(),
        #"statistics" => Vector{String}(),
        "aliases" => ["calculate_weights_climwip"],
        #"timeranges" => Vector{String}()
    )
);




# 1. Compute Weights
# configure parameters for computing weights
config_weights = sw.ConfigWeights(
    performance = Dict("tas_CLIM"=>1, "pr_CLIM"=>2, "psl_CLIM"=>1),
    independence = Dict("tas_CLIM"=>0.5, "pr_CLIM"=>0.25, "psl_CLIM"=>0),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    ref_period = "1995-2014"
);


weights = sw.getOverallWeights(data_weights, config_weights);
# make some Plots
sw.plotWeights(weights.overall)

wP_by_var = dropdims(reduce(+, weights.performance_all, dims=:diagnostic), dims=:diagnostic)
figs_performance = sw.plotPerformanceWeights(wP_by_var)

figs_Di = sw.plotPerformanceWeights(wP_by_var; wP_combined=weights.performance, isBarPlot=false);
f_Di = sw.plotPerformanceWeights(wP_by_var; isBarPlot=false)

figs_Sij = sw.plotIndependenceWeights(weights.independence)
wI_by_var = dropdims(reduce(+, weights.independence_all, dims=:diagnostic), dims=:diagnostic)
figs_wI = sw.plotIndependenceWeights(wI_by_var)



# 2. apply computed weights - Temperature map plots
data_temp_map_future = sw.loadData(
    config_path,
    base_path, 
    dir_per_var=false,
    constraints = sw.DataConstraint(
        tasks = ["weighted_temperature_map_future"]
    )
);
data_temp_map_reference = sw.loadData(
    config_path,
    base_path, 
    dir_per_var=false,
    constraints = sw.DataConstraint(
        tasks = ["weighted_temperature_map_reference"]
    )
);

# compute weighted averages and plot results
# make sure that same models for reference and time period of interest are used
sw.keepSharedModelData!(data_temp_map_future.models, data_temp_map_reference.models)

data_ref = data_temp_map_reference.models["historical-rcp85_CLIM_tas_1995-2014#weighted_temperature_map_reference"]
data_future = data_temp_map_future.models["historical-rcp85_CLIM_tas_2081-2100#weighted_temperature_map_future"]
# to align with original data
data_ref = data_ref[lat = Where(x -> x <= 68.75)];
data_future = data_future[lat = Where(x -> x <= 68.75)];

weights_all_members = sw.makeWeightPerEnsembleMember(weights.overall);
weighted_ref = sw.computeWeightedAvg(data_ref, weights_all_members)
weighted_future = sw.computeWeightedAvg(data_future, weights_all_members)
diff_ww = weighted_future .- weighted_ref;
diff_ww = sw.sortLongitudesWest2East(diff_ww)

f1 = sw.plotMeansOnMap(diff_ww, "Weighted mean temp. change 2081-2100 - 1995-2014", ColorSchemes.Reds.colors)



function compareToOrigData(data, data_orig)
    mat = isapprox.(data, data_orig, atol=10^-4);
    approxEqual = all(x -> ismissing(x) || x, mat)
    @assert approxEqual
end

begin
    # compare to original data (must have adapted the latitudes above)
    # TODO: check why in my data more missing values!
    data_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_weighted_map.nc");
    data_ww_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(diff_ww, data_ww_orig)
    d = DimArray(data_ww_orig, (Dim{:lon}(Array(dims(diff_ww, :lon))), Dim{:lat}(Array(dims(diff_ww, :lat)))))
    sw.plotMeansOnMap(d, "Weighted mean temp. change 2081-2100 - 1995-2014", ColorSchemes.Reds.colors);
end

unweighted_ref = sw.computeWeightedAvg(data_ref);
unweighted_future = sw.computeWeightedAvg(data_future);
diff_uu = unweighted_future .- unweighted_ref;
diff_uu = sw.sortLongitudesWest2East(diff_uu);
diff = diff_ww .- diff_uu

f2 = sw.plotMeansOnMap(
    diff,
    "Weighted minus unweighted mean temp. change: 2081-2100 minus 1995-2014",
    ColorSchemes.Reds.colors
)

# compare to orig data (TODO: again check missing values!)
begin
    data_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_difference_map.nc");
    data_wu_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    compareToOrigData(diff, data_wu_orig)
end


# 3. Apply computed weights - Temperature graph plots
data_temp_graph = sw.loadData(
    config_path,
    base_path, 
    dir_per_var=false,
    constraints = sw.DataConstraint(
        tasks = ["weighted_temperature_graph"]
    )
);
data_graph = data_temp_graph.models["historical-rcp85_ANOM_tas_1960-2100#weighted_temperature_graph"];
#data_graph = sw.averageEnsembleVector(data_graph, true);

uncertainties = sw.getUncertaintyRanges(data_graph, weights_all_members);
weighted_avg = sw.computeWeightedAvg(data_graph, weights_all_members);
unweighted_avg = sw.computeWeightedAvg(data_graph);

f3 = sw.plotTempGraph(
    data_graph, 
    (weighted=weighted_avg, unweighted=unweighted_avg),
    uncertainties, 
    "Temperature anomaly relative to 1981-2010", 
    ""
)
tas_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_graph/weighted_temperature_graph/temperature_anomalies.nc")["tas"];
@assert tas_orig == data_graph


# TODO: check differences in data (missing for BNU-ESM?!)
# bnu_future_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/preproc/weighted_temperature_map/tas_CLIM_future/CMIP5_BNU-ESM_Amon_historical-rcp85_r1i1p1_tas_2081-2100.nc")["tas"];
# bnu_ref_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/preproc/weighted_temperature_map/tas_CLIM_reference/CMIP5_BNU-ESM_Amon_historical-rcp85_r1i1p1_tas_1995-2014.nc")["tas"];
# bnu_future_orig[:,:] .=== bnu_future.data
# bnu_ref_orig .=== bnu_ref.data


# TODO: Plot ensemble spread for some models with >1 ensemble member
# data = modelDataRef["CLIM"]["tas"]
# models_kept = ["EC-Earth3", "ACCESS-CM2"]
# indices = findall(m -> m in models_kept, Array(dims(data, :model)))
# shared_models = data.metadata["full_model_names"][indices]
# data = sw.getModelSubset(data, shared_models)
# sw.plotEnsembleSpread(data, 7.5, 82.5)
