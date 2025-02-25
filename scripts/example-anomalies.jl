import ModelWeights as mw
using NCDatasets
using DimensionalData
using Setfield
using CairoMakie
using ColorSchemes
using Colors


# LGM simulations
path_config = "/albedo/home/brgrus001/ModelWeights/configs/examples/example-anomalies-lgm-piControl.yml";
lgm_data =  mw.loadDataFromYAML(path_config, subset=Dict("level_shared_models" => mw.MEMBER));
mw.computeAnomalies!(lgm_data, "tas_CLIM_lgm", "tas_CLIM_piControl")

colorrange = reverse(Colors.colormap("RdBu", logscale=false, mid=0.25));
f1 = mw.makeSubplots(lgm_data.map["tas_ANOM_lgm"].data, (nrows=3, ncols=5);
    figsize = (800,800).*(5,3), color_range_limits = (-30, 10), title="LGM minus piControl",
    colors = colorrange[2:end-1], low_clip = colorrange[1], high_clip = colorrange[end] 
)
save("plots/anomalies/lgm_anomalies.png", f1)

# Historical simulations
path_config = "/albedo/home/brgrus001/ModelWeights/configs/examples/example-anomalies-historical.yml";
historical_data =  mw.loadDataFromYAML(path_config, subset=Dict("level_shared_models" => mw.MEMBER));
# take only the exact same models for all variables 
df_historical = mw.subsetModelData(historical_data)
mw.computeAnomalies!(df_historical, "tas_CLIM_historical3", "tas_CLIM_historical0");
mw.computeAnomalies!(df_historical, "tos_CLIM_historical3", "tos_CLIM_historical0");


# make plots for all members of a single model
members = Array(dims(df_historical.map["tas_ANOM_historical3"].data, :member));
indices = findall(x -> startswith(x, "ACCESS-CM2#"), members);
data = df_historical.map["tas_ANOM_historical3"].data[member = Where(x -> x in members[indices])]
colorrange = reverse(colormap("RdBu", logscale=false));
f2 = mw.makeSubplots(data, (nrows=2, ncols=5); figsize= (800,600).*(5,3), 
    color_range_limits = (-2.5, 2.5), 
    title="Mean Near-Surface Air Temperature Anomalies 1991-2014 wrt 1850-1900",
    colors = colorrange[2:end-1], low_clip = colorrange[1], high_clip = colorrange[end] 
    )
save("plots/anomalies/historical_anomalies.png", f2)

# add land/sea masks
mw.addMasks!(df_historical, "orog_none_historical")


# land/sea mask
ocean_mask = mw.getOceanMask(df_historical.map["orog_none_historical"].data);
land_mask = ocean_mask .== false;

# plot masks for single model (not identical for all due to different missing values/possibly different orog values)
f_ocean = Figure();
mw.plotMeansOnMap!(f_ocean, ocean_mask[:,:,1], "ocean mask")
f_land = Figure();
mw.plotMeansOnMap!(f_land, land_mask[:,:,1], "land mask")


# compute global mean temperature anomaly for every member for land/sea
tas_anom_data = df_historical.map["tas_ANOM_historical3"].data
# sea
ocean_tas = similar(tas_anom_data, Union{Missing, Float32});
ocean_tas .= tas_anom_data;
ocean_tas[land_mask] .= missing;

# plot ocean data for a single model
model = "TaiESM1#r1i1p1f1_gn"
f1 = Figure();
cmap = reverse(Colors.colormap("RdBu", mid=2/3));
mw.plotMeansOnMap!(f1, ocean_tas[member = At(model)], "Ocean Anomalies 1991-2014 minus 1850-1900: $(model)";
    colors = cmap[2:end-1],
    color_range = (-1, 2), 
    high_clip = cmap[end], low_clip = cmap[1]
)
f1

# land
land_tas = similar(tas_anom_data, Union{Missing, Float32});
land_tas .= tas_anom_data;
land_tas[ocean_mask] .= missing;
# plot land data for a single model
f2= Figure();
mw.plotMeansOnMap!(f2, land_tas[member = At(model)], "Land Anomalies 1991-2014 minus 1850-1900: $(model)";
    colors = cmap[2:end-1],
    color_range = (-1, 2), 
    high_clip = cmap[end], low_clip = cmap[1]
)
f2



gms_tas_everywhere = mw.getGlobalMeans(tas_anom_data);
gms_tas_ocean = mw.getGlobalMeans(ocean_tas);
gms_tas_land = mw.getGlobalMeans(land_tas);
gms_ratio = gms_tas_land ./ gms_tas_ocean;

non_missing_indices = findall(x -> !ismissing(x), gms_ratio)

f3, ax = mw.makeScatterPlot(
    gms_tas_everywhere[non_missing_indices], gms_ratio[non_missing_indices]; 
    captions = (
        x="mean temperature anomaly all locations", 
        y="ratio of global mean temperature anomaly on land vs. ocean",
        title="Historical model simulations for period 1991-2014 minus 1850-1900"),
    add_lines = false
)
f3
save("plots/anomalies/gm_anomalies_land_vs_ocean.png", f3)

# plot tas vs. ratio of tas on land/tos (proxy for ratio)
tos_anom_data = df_historical.map["tos_ANOM_historical3"].data;
gms_tos = mw.getGlobalMeans(tos_anom_data);

f4, ax = mw.makeScatterPlot(
    gms_tas_land, gms_tos; 
    captions = (
        x="global mean anomaly tas on land", 
        y="global mean anomaly tos (sst, only defined on sea)",
        title="Historical model simulations for period 1991-2014 minus 1850-1900"),
    add_lines = false
)
f4
save("plots/anomalies/gm_anomalies_tas-land_vs_tos-sea.png", f4)
