using NetCDF
using PyCall, PyPlot
using Dates

@pyimport mpl_toolkits.basemap as basemap
@pyimport numpy as np
@pyimport textwrap

include(joinpath(@__DIR__, "..", "..", "src", "load-data-utils.jl"));
include(joinpath(@__DIR__, "..", "..", "src", "plot-utils.jl"));

workDir = joinpath(PATH_TO_WORK_DIR, "weighted_temperature_map", "weighted_temperature_map");
targetDir = joinpath("reproduce-climwip-figs", "plots-replicated-with-julia", "weighted_temperature_map_pycall");
mkpath(targetDir)



function plotWeightedTemperatureMap(workDir, filenameData, textTitle, vmin=-1, vmax=1, cmap="RdBu_r")
    data = NetCDF.ncread(joinpath(workDir, filenameData), "__xarray_dataarray_variable__");
    latY = NetCDF.ncread(joinpath(workDir, filenameData), "lat");
    lonX = NetCDF.ncread(joinpath(workDir, filenameData), "lon");

    lonTicks = -10:10:40;
    latTicks = 30:10:80;
    lonLabels = longitude2EastWest.(lonTicks);
    latLabels = latitude2NorthSouth.(latTicks);

    # Create the figure and axis
    fig, ax = PyPlot.subplots()

    ax.set_xticks(lonTicks)
    ax.set_yticks(latTicks)
    ax.set_xticklabels(lonLabels)
    ax.set_yticklabels(latLabels)


    # Create a Basemap instance
    # llcrnrlon: longitude of lower left hand corner of the desired map domain (degrees)
    m = basemap.Basemap(projection="cyl", llcrnrlat=latTicks[1], urcrnrlat=last(latTicks), llcrnrlon=lonTicks[1], urcrnrlon=last(lonTicks),
                        resolution="i", ax=ax, suppress_ticks = false)
    m.drawcoastlines()

    # Create a meshgrid for longitude and latitude
    lon, lat = np.meshgrid(lonX, latY)

    # Plot the data using pcolormesh
    cs = m.pcolormesh(lon, lat, data', cmap=cmap, shading="auto", vmin=vmin, vmax=vmax);
    cbar = m.colorbar(cs, location="right", pad="10%")

    # Add title and labels
    wrapped_title = textwrap.fill(textTitle, width=40)
    title(wrapped_title)
    xlabel("Longitude")
    ylabel("Latitude")

    return fig
end


fig1 = plotWeightedTemperatureMap(workDir, "temperature_change_difference_map.nc", "Weighted minus unweighted mean temperature change: 2081-2100 minus 1995-2014 (°C)")
savefig(joinpath(targetDir, getCurrentTime() * "_weighted_temperature_map_pycall.png"));

fig2 = plotWeightedTemperatureMap(workDir, "temperature_change_weighted_map.nc", "Weighted mean temperature change: 2081-2100 minus 1995-2014 (°C)", 2.5, 6.5, "Reds")
savefig(joinpath(targetDir, getCurrentTime() * "_weighted_mean_temp_change_pycall.png"));

