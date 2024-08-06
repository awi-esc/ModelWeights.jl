using NetCDF
using Dates
using CairoMakie
using TextWrap
using GeoMakie

include(joinpath(@__DIR__, "..", "..", "src", "load-data-utils.jl"));
include(joinpath(@__DIR__, "..", "..", "src", "plot-utils.jl"));

workDir = joinpath(PATH_TO_WORK_DIR, "weighted_temperature_map", "weighted_temperature_map");
targetDir = joinpath(pwd(),"reproduce-climwip-figs", "plots-replicated-with-julia", "weighted_temperature_map");
mkpath(targetDir)


function plotWeightedTemperatureMap(workDir, filenameData, textTitle, colors, vmin=-1, vmax=1)
    data = NetCDF.ncread(joinpath(workDir, filenameData), "__xarray_dataarray_variable__");
    latY = NetCDF.ncread(joinpath(workDir, filenameData), "lat");
    lonX = NetCDF.ncread(joinpath(workDir, filenameData), "lon");

    lonTicks = -10:10:40;
    latTicks = 30:10:80;
    lonLabels = longitude2EastWest.(lonTicks);
    latLabels = latitude2NorthSouth.(latTicks);

    # Create the figure and axis
    fig = Figure();
    ax = Axis(fig[1,1], 
        title = wrap(textTitle, width=40),
        xlabel = "Longitude",
        ylabel = "Latitude",
        xticks = (lonTicks, lonLabels),
        yticks = (latTicks, latLabels),
        limits = ((lonTicks[1], last(lonTicks)), (latTicks[1], last(latTicks)))
    );
    hm = heatmap!(ax, lonX, latY, data, colormap = colors, alpha=0.8, colorrange=(vmin, vmax));
    lines!(GeoMakie.coastlines(),color=:black);
    Colorbar(fig[1,2], hm)
    return fig
end


fig1 = plotWeightedTemperatureMap(
    workDir, 
    "temperature_change_difference_map.nc", 
    "Weighted minus unweighted mean temperature change: 2081-2100 minus 1995-2014 (°C)",
    reverse(ColorSchemes.RdBu.colors)
)
save(joinpath(targetDir, getCurrentTime() * "_weighted_temperature_map.png"), fig1);

fig2 = plotWeightedTemperatureMap(
    workDir, 
    "temperature_change_weighted_map.nc",
    "Weighted mean temperature change: 2081-2100 minus 1995-2014 (°C)", 
    ColorSchemes.Reds.colors,
    2.5, 
    6.5
)
save(joinpath(targetDir, getCurrentTime() * "_weighted_mean_temp_change.png"), fig2);

