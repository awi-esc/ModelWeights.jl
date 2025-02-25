using ColorSchemes
using Dates
using TextWrap
using Statistics
using DimensionalData

using CairoMakie  
using GeoMakie  

CairoMakie.activate!(type = "svg")

function getFigure(figsize, fontsize)
    size_pt = 72 .* figsize
    fig=Figure(size= size_pt, fontsize=fontsize)
    return fig
end

function savePlot(fig, target_dir::String, target_fn::String)
    if !isdir(target_dir)
        mkpath(target_dir)
    end
    path_to_target = joinpath(target_dir, target_fn);
    save(path_to_target, fig);
    @info "saved plot to " path_to_target
end

function plotDistMatrices(distMat, diagnostic, models, modelRefs)
    fig = Figure();
    xs = 1:length(models)
    ax = Axis(
        fig[1,1],
        xlabel = "Model",
        ylabel = "Model",
        xticks = (xs, models), 
        yticks = (xs, modelRefs),
        xticklabelrotation = pi/4,
        title = "Distance matrix " * diagnostic,
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
    sortLongitudesWest2East(data::DimArray)

Arrange 'data' such that western latitudes come first, then eastern latitudes.
"""
function sortLongitudesWest2East(data::DimArray)
    lon = dims(data, :lon)
    east = lon[lon .< 180]
    west = lon[lon .>= 180]
    sorted_lon = [Array(west); Array(east)];

    # necessary to specify that lookup dimension values aren't sorted anymore!
    # otherwise Selector At won't work! does seem to work don't know TODO: CHECK THIS!!
    lookup_lon = Lookups.Sampled(
        sorted_lon;
        span=Lookups.Irregular(minimum(lon), maximum(lon)), 
        # order=Lookups.Unordered()
        order = Lookups.ForwardOrdered()
    )
    data = data[lon=At(sorted_lon)]
    data = DimensionalData.Lookups.set(data, lon = lookup_lon)
    return data
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

"""
    kelvinToCelsius(data::DimArray)

Return DimensionalData.DimArray with values that were given in Kelvin converted 
to Degree Celsius.

# Arguments:
- `data`:
"""
function kelvinToCelsius(data::DimArray)
    units =  data.metadata["units"]
    if isa(units, String) && units == "K"
        data = data .- 273.15
    elseif isa(units, Vector)
        indices = findall(units .!= "K")
        if hasdim(data, :member)
            data[member = indices] =  data[member = indices] .- 273.15
        else
            data[model = indices] =  data[model = indices] .- 273.15
        end
    end
    data.metadata["units"] = "degC"
    return data
end


"""
    kelvinToCelsius!(datamap::DataMap)

Modify entries of `datamap` such that all data is given in Degree Celsius (instead) 
of Kelvin.

# Arguments:
- `datamap`:
"""
function kelvinToCelsius!(datamap::DataMap)
    for (id, da) in datamap.map
        df = @set da.data = kelvinToCelsius(da.data)
        addToMap!(datamap, df)
    end
    return nothing
end


"""
    makeSubplots
    
# Arguments:

"""
function makeSubplots(
    data::DimArray, grid::NamedTuple{(:nrows, :ncols), Tuple{Int, Int}}; 
    fontsize=36, figsize=(800,600), title="",
    colors=nothing, color_range_limits=nothing, high_clip=(1,0,0), low_clip=(0,0,1)
)
    models = hasdim(data, :member) ? Array(dims(data, :member)) : 
        (hasdim(data, :model) ? Array(dims(data, :model)) : nothing)
    if isnothing(models)
        throw(ArgumentError("subplots only possible for data with dimension :model or :member"))
    end
    # models = reshape(models, grid...)
    fig = Figure(size = figsize, fontsize=fontsize)
    Label(fig[0, 1:grid.ncols], title, fontsize = 1.5 * fontsize, halign=:center, font=:bold)

    nb_subplots = length(models)

    for idx_plot in 1:nb_subplots
        row = ceil(Int, idx_plot / grid.ncols)
        col_temp = idx_plot % grid.ncols
        col = col_temp==0 ? grid.ncols : col_temp
        pos = (x=row, y=col)
        pos_legend = idx_plot == nb_subplots ? (x=1:row, y=grid.ncols+1) : nothing
        model = models[idx_plot]
        if hasdim(data, :member)
            plotMeansOnMap!(
                fig, data[member = At(model)], "$model";
                colors=colors, high_clip=high_clip, low_clip=low_clip, 
                color_range=color_range_limits, pos=pos, pos_legend=pos_legend)
        else
            plotMeansOnMap!(
                fig, data[model = At(model)], "$model";
                colors=colors, high_clip=high_clip, low_clip=low_clip,
                color_range=color_range_limits, pos=pos, pos_legend=pos_legend)
        end
    end
    return fig
end