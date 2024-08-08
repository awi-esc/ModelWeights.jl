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


function plotDistMatrices(distMat, climateVar, models, modelRefs)
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



""" longitude2EastWest(lon::Number)

Convert longitudes from -180° to 180° into 0° to 180° East/West.
"""
function longitude2EastWest(lon::Number)
    return lon > 0 ? "$(lon)°E" : "$(abs(lon))°W"
end

""" latitude2NorthSouth(lat::Number)

Convert latitudes from -90° to 90° into 0° to 90° North/South.
"""
function latitude2NorthSouth(lat::Number)
    return lat < 0 ? "$(abs(lat))°S" : "$(lat)°N"
end

""" 
    lon360to180(lon::Number)
    
Convert longitudes measured from 0° to 360° into  -180° to 180° scale.
"""
function lon360to180(lon::Number)
    return lon > 179 ? lon-360 : lon
end

""" 
    lon360to180(lon::Number)

Convert longitudes measured from -180° to 180° into 0° to 360° scale.
"""
function lon180to360(lon::Number)
    return ifelse(lon < 0, lon + 360, lon)
end


"""
    convertKgsToSv!(vec:DimArray)

Convert data given in unit 'kg s-1' into Sverdrups (Sv).
"""
function convertKgsToSv!(data::DimArray)
    if data.metadata["units"] != "kg s-1"
        throw(ArgumentError("The unit of the data should be 'kg s-1', but it is " * data.metadata["units"]))
    end
    data[1:end] = data .* (10^-9);
    data.metadata["units"] = "Sv";
    return nothing
end



""" plotMeansOnMap(means::DimArray, title::String)
    
Plot contours of world with an overlayed heatmap that shows the data which correspond to
mean value for each position in considered grid. 
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

Find the grid point in given `longitudes`, `latitudes` - grid closest to `location`.

# Arguments
- `location::Dict`: 'lon', 'lat' of position for which closest grid point is returned, `lon`
must be given from -180° to 180°
- `longitudes::Vector`: grid longitudes measured from -180° to 180°
- `latitudes::Vector`: grid latitudes measured from -90° to 90°
"""
function getClosestGridPoint(location::Dict, longitudes::Vector, latitudes::Vector)
    idx_lat = argmin(abs.(latitudes .- location["lat"]));
    idx_lon = argmin(abs.(longitudes .- location["lon"]));
    lat = latitudes[idx_lat];
    lon = longitudes[idx_lon];
    return Dict([("name", location["name"]), ("lon", lon), ("lat", lat)])
end

""" plotHistAtPos(data::DimArray, location::Dict)

Plot histogram of all data for a specific `location`.

# Arguments:
- `data::DimArray`: dimensions 'lon' (from -180° to 180°), 'lat' (-90° to 90°)
- `location::Dict`: keys 'name', 'lon', 'lat'
"""
function plotHistAtPos(data::DimArray, location::Dict)
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
    ax = Axis(fig[1,1], title = join([t1, t2], "\n"), xlabel = data.metadata["units"])
    hist!(ax, vec(data_loc))
    return fig
end


""" plotAMOC(data::DimArray)

Plot AMOC strength for variable "amoc".
"""
function plotAMOC(data::DimArray)
    fig = getFigure((8, 5), 12);
    t = "variable: amoc, experiment: " * data.metadata["experiment_id"];

    unit = data.metadata["units"];
    ax = Axis(fig[1, 1], title=join(["AMOC strength (at 26.5°N)", t], "\n"),  xlabel = "Transport in " * unit)
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

