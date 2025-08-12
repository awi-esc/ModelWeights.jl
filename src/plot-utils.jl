function getFigure(figsize, fontsize)
    size_pt = 72 .* figsize
    fig = Figure(size = size_pt, fontsize = fontsize)
    return fig
end

function savePlot(fig, target_path::String; overwrite::Bool=false)
    target_path = overwrite ? target_path : Data.individuatePath(target_path)
    save(target_path, fig)
    @info "saved plot to " target_path
end

function savePlot(fig, target_dir::String, target_fn::String; overwrite::Bool=false)
    savePlot(fig, joinpath(target_dir, target_fn); overwrite)
end

function plotDistMatrices(distMat, diagnostic, models, modelRefs)
    fig = Figure()
    xs = 1:length(models)
    ax = Axis(
        fig[1, 1],
        xlabel = "Model",
        ylabel = "Model",
        xticks = (xs, models),
        yticks = (xs, modelRefs),
        xticklabelrotation = pi / 4,
        title = "Distance matrix " * diagnostic,
        yreversed = true,
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
    convertKgsToSv!(vec:YAXArray)

Convert data given in unit 'kg s-1' into Sverdrups (Sv).
"""
function convertKgsToSv!(data::YAXArray)
    if data.properties["units"] != "kg s-1"
        throw(ArgumentError("Required unit: 'kg s-1', found: $(data.properties["units"])"))
    end
    data[1:end] = data .* (10^-9)
    data.properties["units"] = "Sv"
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
    idx_lat = argmin(abs.(latitudes .- location["lat"]))
    idx_lon = argmin(abs.(longitudes .- location["lon"]))
    lat = latitudes[idx_lat]
    lon = longitudes[idx_lon]
    return Dict([("name", location["name"]), ("lon", lon), ("lat", lat)])
end


"""
    makeSubplots
    
# Arguments:

"""
function makeSubplots(
    data::AbstractArray,
    grid::NamedTuple{(:nrows, :ncols),Tuple{Int,Int}};
    fontsize = 12,
    figsize = (600, 450),
    title = "",
    colors = nothing,
    color_range_limits = nothing,
    high_clip = (1, 0, 0),
    low_clip = (0, 0, 1),
    xlabel = "Longitude",
    ylabel = "Latitude",
    xlabel_rotate = pi / 4,
)
    models =
        hasdim(data, :member) ? Array(dims(data, :member)) :
        (hasdim(data, :model) ? Array(dims(data, :model)) : nothing)
    if isnothing(models)
        throw(
            ArgumentError(
                "subplots only possible for data with dimension :model or :member",
            ),
        )
    end
    # models = reshape(models, grid...)
    fig = Figure(size = figsize, fontsize = fontsize)
    Label(
        fig[0, 1:grid.ncols],
        title,
        fontsize = 1.5 * fontsize,
        halign = :center,
        font = :bold,
    )

    nb_subplots = length(models)

    for idx_plot = 1:nb_subplots
        row = ceil(Int, idx_plot / grid.ncols)
        col_temp = idx_plot % grid.ncols
        col = col_temp == 0 ? grid.ncols : col_temp
        pos = (x = row, y = col)
        pos_legend = idx_plot == nb_subplots ? (x = 1:row, y = grid.ncols + 1) : nothing
        model = models[idx_plot]
        if hasdim(data, :member)
            plotValsOnMap!(
                fig,
                data[member=At(model)],
                "$model";
                colors = colors,
                high_clip = high_clip,
                low_clip = low_clip,
                color_range = color_range_limits,
                pos = pos,
                pos_legend = pos_legend,
                xlabel = xlabel,
                ylabel = ylabel,
                xlabel_rotate = xlabel_rotate,
            )
        else
            plotValsOnMap!(
                fig,
                data[model=At(model)],
                "$model";
                colors = colors,
                high_clip = high_clip,
                low_clip = low_clip,
                color_range = color_range_limits,
                pos = pos,
                pos_legend = pos_legend,
                xlabel = xlabel,
                ylabel = ylabel,
                xlabel_rotate = xlabel_rotate,
            )
        end
    end
    return fig
end
