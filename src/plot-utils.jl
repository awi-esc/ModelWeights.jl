using ColorSchemes
using Dates
using TextWrap
using Statistics
using DimensionalData

using CairoMakie  
using GeoMakie  
CairoMakie.activate!(type = "svg")

function getFigure(figsize, fontsize)
    size_inches = figsize#(8, 5)
    size_pt = 72 .* size_inches
    fig=Figure(size= size_pt, fontsize=fontsize)
    return fig
end

"""
"""
function plotDistMatrices(distMat, climateVar, models, modelRefs)
    # size_inches = (8, 5)
    # size_pt = 72 .* size_inches
    # fig=Figure(size= size_pt, fontsize=12)
    fig = Figure();
    xs = 1:length(models)
    ax = Axis(
        fig[1,1],
        xticks = (xs, models), 
        yticks = (xs, modelRefs),
        xticklabelrotation = pi/4,
        title = "Distance matrix for " * climateVar,
        yreversed = true
        )
    hm = heatmap!(ax, distMat', colormap = ColorSchemes.YlGn_4.colors)
    Colorbar(fig[1, 2], hm)
    return fig
end



""" longitude2EastWest(lon)
    converts longitudes from -180° to 180° into 0° to 180° East/West
"""
function longitude2EastWest(lon)
    return lon > 0 ? "$(lon)°E" : "$(abs(lon))°W"
end

""" latitude2NorthSouth(lat)
    converts latitudes from -90° to 90° into 0° to 90° North/South.
"""
function latitude2NorthSouth(lat)
    return lat < 0 ? "$(abs(lat))°S" : "$(lat)°N"
end

""" 
    lon360to180(lon::Number)
    
    converts longitudes measured from 0° to 360° into degrees measured from -180° to 180°.
"""
function lon360to180(lon::Number)
    return lon > 179 ? lon-360 : lon
end

""" 
    lon360to180(lon::Number)

    converts longitudes measured from -180° to 180° into degrees measured from 0° to 360°.
"""
function lon180to360(lon::Number)
    return ifelse(lon < 0, lon + 360, lon)
end




""" plotMeansOnMap(means::DimArray, title::String)
    plots contours of world with an overlayed heatmap that shows the data which correspond to
    mean value for each position in considered grid. 
        
    Longitudes are supposed to be given as measured from -180° to 180°.
"""
function plotMeansOnMap(means::DimArray, title::String)
    dims_lat = Array(dims(means, :lat));
    dims_lon = Array(dims(means, :lon));
    if any(x -> x > 179, dims_lon)
        dims_lon = lon360to180.(dims_lon);
    end
    
    # scaling plot 
    lon_min, lon_max = -180, 180;
    lat_min, lat_max = -90, 90;
    lon = range(lon_min, stop=lon_max, length=length(dims_lon));
    lat = range(lat_min, stop=lat_max, length=length(dims_lat));

    # axis ticks
    lonLabels = SimilarityWeights.longitude2EastWest.(dims_lon);
    latLabels = SimilarityWeights.latitude2NorthSouth.(dims_lat);
    # just use roughly 10 ticks
    step_lon = Int(round(length(lonLabels)/10));
    step_lat = Int(round(length(latLabels)/10));

    fig = Figure();
    ax = Axis(fig[1,1], 
        title = TextWrap.wrap(title, width=40),
        xlabel = "Longitude",
        ylabel = "Latitude",
        xticklabelrotation = pi/4,
        xticks = (dims_lon[1 : step_lon : end], lonLabels[1 : step_lon : end]),
        yticks = (dims_lat[1 : step_lat : end], latLabels[1 : step_lat : end]),
        limits = ((lon_min, lon_max), (lat_min, lat_max))
        );
    lines!(GeoMakie.coastlines(); color=:black);
    hm = heatmap!(ax, lon, lat, Array(means), alpha=0.8);
    Colorbar(fig[1,2], hm);
    return fig
end



""" getClosestGridPoint(location::Dict, longitudes::Vector, latitudes::Vector)

    location is a dictionary with keys: lon, lat mapping to the location for which the
    closest point on the given grid is returned.

Longitudes are assumed to be measured from -180 to 180 and latitudes from -90 to 90 degrees.
Make sure that location refers to this scale, too. E.g. 96°W has to be given as -96 longitude.
"""
function getClosestGridPoint(location::Dict, longitudes::Vector, latitudes::Vector)
    idx_lat = argmin(abs.(latitudes .- location["lat"]));
    idx_lon = argmin(abs.(longitudes .- location["lon"]));
    lat = latitudes[idx_lat];
    lon = longitudes[idx_lon];
    return Dict([("name", location["name"]), ("lon", lon), ("lat", lat)])
end

""" plotHistAtPos(data::DimArray, location::Dict, unit::String)

    plots histogram of all model data for a specific location.

data must have dimnensions 'lon', 'lat', which are given from -180 to 180 (lon) and from -90 to 90 (lat).
location is a Dict with keys 'name', 'lon', 'lat'
(unit should be retrieved via metadata from data in the future)
"""
function plotHistAtPos(data::DimArray, location::Dict, unit::String="")
    coords = getClosestGridPoint(location, Array(dims(data, :lon)), Array(dims(data, :lat)));
    # no idea why At doesnt work for negative values at longitude dimension...
    # At would drop dimensions directly..
    data_loc = data[lat=Where(elem-> elem == coords["lat"]), lon=Where(elem -> elem == coords["lon"])];
    data_loc = dropdims(data_loc, dims=:lon)
    data_loc = dropdims(data_loc, dims=:lat)

    grid_lat = latitude2NorthSouth(coords["lat"])
    grid_lon = longitude2EastWest(coords["lon"])
    t1 = "Variable: " * data.metadata["variable_id"] * " Experiment: " * data.metadata["experiment_id"];
    t2 = "near " * location["name"] * "(" * grid_lat * "," * grid_lon * ")";
    
    fig = getFigure((14,10), 16);
    ax = Axis(fig[1,1], title = join([t1, t2], "\n"), xlabel = unit)
    hist!(ax, vec(data_loc))
    return fig
end


""" plotAMOC(data)
    plots AMOC strength for variable "amoc" derived by ESMValTool recipe

unit in ylabel should be done via metadata
"""
function plotAMOC(data)
    fig = getFigure((8, 5), 12);
    t = "variable: amoc, experiment: " * data.metadata["experiment_id"];

    ax = Axis(fig[1, 1], title=join(["AMOC strength (at 26.5°N)", t], "\n"),  ylabel = "Transport in kg s^-1")
    if length(dims(data)) == 1
        # data for full period
        hist!(ax, Array(data))
    else
        # seasonal data
        for i in 1:4
            hist!(ax, Array(data[i,:]), scale_to=-0.6, offset=i, direction=:x)
        end
        ax.xlabel = "season"
    end
    return fig
end

