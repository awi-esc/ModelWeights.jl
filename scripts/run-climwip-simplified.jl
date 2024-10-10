import SimilarityWeights as sw
using NCDatasets
using ColorSchemes
using CairoMakie
using GeoMakie
using TextWrap
using DimensionalData

begin
    path_config_weights = "configs/climwip_simplified_weights.yml";
    config_weights = sw.validateConfig(path_config_weights);

    # get overall weights
    weights = sw.getOverallWeights(config_weights);
    sw.plotWeights(weights)

    ##### Performance and Independence weights #####
    modelDataFull, modelDataRef = sw.loadModelData(config_weights);
    obsData = sw.loadDataFromConfig(
        config_weights, 
        config_weights.name_ref_period, 
        config_weights.obs_data_name
    );    
    wP = sw.getPerformanceWeights(
        modelDataRef["CLIM"], 
        obsData["CLIM"], 
        config_weights.weights_variables["performance"]
    )
    Di = sw.reduceGeneralizedDistancesVars(wP);
    # or set summarize_variables to true: 
    # Di = sw.getPerformanceWeights(
    #     modelDataRef["CLIM"], 
    #     obsData["CLIM"], 
    #     config_weights.weights_variables["performance"], 
    #     true
    # )
    wI = sw.getIndependenceWeights(
        modelDataRef["CLIM"], 
        config_weights.weights_variables["independence"]
    );
    Sij = sw.reduceGeneralizedDistancesVars(wI);
    # or for computing Sij, set summarize_variables=true

    figs_wP = sw.plotPerformanceWeights(wP)
    figs_Di = sw.plotPerformanceWeights(wP, Di, false)
    figs_Di = sw.plotPerformanceWeights(wP, nothing, false)

    figs_wI = sw.plotIndependenceWeights(wI)
    figs_Sij = sw.plotIndependenceWeights(Sij)
end

##### Temperature graph #####
begin
    weights = sw.loadWeightsAsDimArray("/albedo/work/projects/p_forclima/britta/similarityweights-output/climwip-simplified/2024-10-09_08_25/weights.nc");
    weights_all_members = sw.makeWeightPerEnsembleMember(weights);
    
    path_config_temp_graph = "configs/climwip_simplified_temperature_graph.yml";
    config_graph = sw.validateConfig(path_config_temp_graph);

    modelDataFull, _ = sw.loadModelData(config_graph);

    tas_data_all_members = modelDataFull["ANOM"]["tas"];
    tas_data = sw.averageEnsembleVector(tas_data_all_members, true);

    # check that data is identical to that from original climwip
    tas_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_graph/weighted_temperature_graph/temperature_anomalies.nc")["tas"];
    @assert tas_orig == tas_data_all_members

    uncertainties = sw.getUncertaintyRanges(tas_data_all_members, weights_all_members);
    weighted_avg = sw.computeWeightedAvg(tas_data_all_members, weights_all_members)
    unweighted_avg = sw.computeWeightedAvg(tas_data_all_members)

    tit = "Temperature anomaly relative to "; #* name_ref_period
    sw.plotTempGraph(
        tas_data_all_members, 
        (weighted=weighted_avg, unweighted=unweighted_avg),
        uncertainties, 
        tit, 
        ""
    )
end

##### Temperature map #####
begin
    weights = sw.loadWeightsAsDimArray("/albedo/work/projects/p_forclima/britta/similarityweights-output/climwip-simplified/2024-10-09_08_25/weights.nc");
    weights_all_members = sw.makeWeightPerEnsembleMember(weights);
        
    path_config_temp_map = "configs/climwip_simplified_temperature_map.yml";
    config_map = sw.validateConfig(path_config_temp_map);
    modelDataFull, modelDataRef = sw.loadModelData(config_map);
    # use same data as in climwip (they only have 20x16 instead of 20x19)
    data_full = modelDataFull["CLIM"]["tas"][lat = Where(x -> x <= 68.75)];
    data_ref =  modelDataRef["CLIM"]["tas"][lat = Where(x -> x <= 68.75)];

    # Compute weighted mean temperature change: 2081-2100 minus 1995-2014 (°C)"
    avgs_full = sw.computeWeightedAvg(data_full, weights_all_members);
    avgs_ref = sw.computeWeightedAvg(data_ref, weights_all_members);
    diff_ww = avgs_full .- avgs_ref;
    diff_ww = sw.sortLongitudesWest2East(diff_ww)
    
    # TODO: diff_ww has more missing data than in original data from climwip
    sw.plotMeansOnMap(diff_ww, "Weighted mean temp. change 2081-2100 - 1995-2014", ColorSchemes.Reds.colors)
    
    # compare to original data (must have adapted the latitudes above)
    data_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_weighted_map.nc");
    data_ww_orig = data_orig["__xarray_dataarray_variable__"][:,:];
    
    mat = isapprox.(diff_ww.data, data_ww_orig, atol=10^-4);
    approxEqual = all(x -> ismissing(x) || x, mat)
    @assert approxEqual

    d = DimArray(data_ww_orig, (Dim{:lon}(Array(dims(diff_ww, :lon))), Dim{:lat}(Array(dims(diff_ww, :lat)))))
    sw.plotMeansOnMap(d, "Weighted mean temp. change 2081-2100 - 1995-2014", ColorSchemes.Reds.colors);
end

# TODO: go through these:
# data_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/work/weighted_temperature_map/weighted_temperature_map/temperature_change_difference_map.nc")
# data_wu_orig = data_orig["__xarray_dataarray_variable__"][:,:]
# diff = modelDataFull["CLIM"]["tas"] .- modelDataRef["CLIM"]["tas"];
# diff_avg = sw.computeWeightedAvg(diff, weights_all_members)
# diff_avg = sw.computeWeightedAvg(sw.averageEnsembleMembers(diff, true), weights)
# df::Dict{String, Dict{String, DimArray}} = Dict("weighted_mean_temp_change" => Dict("tas" => DimArray(diff_ww)))
# sw.plotMeanData(config_map, df)
# averages_full_period = sw.getWeightedAverages(config_map, weights_all_members, "full");
# averages_ref_period = sw.getWeightedAverages(config_map, weights_all_members, "ref");
# diff_ww = averages_full_period["weighted"]["tas"] .- averages_ref_period["weighted"]["tas"];
# diff_wu = averages_full_period["weighted"]["tas"] .- averages_ref_period["unweighted"]["tas"];
# tit = "Weighted minus unweighted mean temperature change: 2081-2100 minus 1995-2014 (°C)";





# colors = reverse(ColorSchemes.RdBu.colors);
# begin
#     #df::Dict{String, Dict{String, DimArray}} = Dict("weighted_mean_temp_change" => Dict("tas" => DimArray(diff)))
#     #sw.plotMeanData(config_map, df)
    
#     plotWeightedTemperatureMap(tit, diff, colors, -1, 1)


#     df2::Dict{String, Dict{String, DimArray}} = Dict("weighted_minus_unweighted_mean_temp_change" => Dict("tas" => DimArray(diff)))
#     # sw.plotMeanData(config_map, df2)
# end

# diff = data_ww_orig
# vmin = -1 
# vmax = 1
# lonTicks = -10:10:40;
# latTicks = 30:10:80;
# lonLabels = sw.longitude2EastWest.(lonTicks);
# latLabels = sw.latitude2NorthSouth.(latTicks);

# lonX = Array(dims(diff, :lon));
# latY = Array(dims(diff, :lat));

# # Create the figure and axis
# fig = Figure();
# ax = Axis(fig[1,1], 
#     title = wrap(tit, width=40),
#     xlabel = "Longitude",
#     ylabel = "Latitude",
#     xticks = (lonTicks, lonLabels),
#     yticks = (latTicks, latLabels),
#     limits = ((lonTicks[1], last(lonTicks)), (latTicks[1], last(latTicks)))
# );
# hm = heatmap!(ax, lonX, latY, parent(diff2), colormap = colors, alpha=0.8, colorrange=(vmin, vmax));
# lines!(GeoMakie.coastlines(),color=:black);
# Colorbar(fig[1,2], hm)
# fig



# check differences in data (missing for BNU-ESM)
bnu_future_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/preproc/weighted_temperature_map/tas_CLIM_future/CMIP5_BNU-ESM_Amon_historical-rcp85_r1i1p1_tas_2081-2100.nc")["tas"];
bnu_ref_orig = NCDataset("/albedo/home/brgrus001/SimilarityWeights/reproduce-climwip-figs/recipe_climwip_test_basic_data/preproc/weighted_temperature_map/tas_CLIM_reference/CMIP5_BNU-ESM_Amon_historical-rcp85_r1i1p1_tas_1995-2014.nc")["tas"];

path_config_temp_map = "configs/climwip_simplified_temperature_map.yml";
config_map = sw.validateConfig(path_config_temp_map);
modelDataFull, modelDataRef = sw.loadModelData(config_map);

bnu_future = modelDataFull["CLIM"]["tas"][model = At("BNU-ESM#r1i1p1")]
bnu_ref = modelDataRef["CLIM"]["tas"][model = At("BNU-ESM#r1i1p1")]


bnu_future_orig[:,:] .=== bnu_future.data
bnu_ref_orig .=== bnu_ref.data