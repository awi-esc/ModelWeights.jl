using ColorSchemes
using Dates
using TextWrap
using Statistics
using DimensionalData

using CairoMakie  
using GeoMakie  
CairoMakie.activate!(type = "svg")


@kwdef struct Target
    directory::String
    filename::String
    save::Bool
end

function getFigure(figsize, fontsize)
    size_pt = 72 .* figsize
    fig=Figure(size= size_pt, fontsize=fontsize)
    return fig
end

function savePlot(fig, target::Union{Target, Nothing}=nothing)
    if !isnothing(target) && target.save
        target_dir = joinpath(target.directory, getCurrentTime());
        if !isdir(target_dir)
            mkpath(target_dir)
        end
        save(joinpath(target_dir, target.filename), fig);
    end
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

function getCurrentTime()
    currentDay = string(Dates.today()) * '_';
    currentTime = Dates.format(Dates.now(), "HH_MM");
    return currentDay * currentTime
end


"""
    convertKgsToSv!(vec:DimArray)

Convert data given in unit 'kg s-1' into Sverdrups (Sv).
"""
function convertKgsToSv!(data::DimArray)
    if data.metadata["units"] != "kg s-1"
        msg = "The unit of the data should be 'kg s-1', but it is " * data.metadata["units"];
        throw(ArgumentError(msg))
    end
    data[1:end] = data .* (10^-9);
    data.metadata["units"] = "Sv";
    return nothing
end


""" getClosestGridPoint(location::Dict, longitudes::Vector, latitudes::Vector)

Find the grid point in grid defined by `longitudes` and `latitudes` that is
closest to `location`.

# Arguments
- `location::Dict`: 'lon', 'lat' of position for which closest grid point is
returned, `lon` must be given from -180° to 180°
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