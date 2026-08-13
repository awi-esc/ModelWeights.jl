function getFigure(figsize, fontsize)
    size_pt = 72 .* figsize
    fig = Figure(size = size_pt, fontsize = fontsize)
    return fig
end

function savePlot(fig, target_path::String; overwrite::Bool=false, verbose::Bool=true)
    target_path = overwrite ? target_path : Data.individuatePath(target_path)
    if !isdir(dirname(target_path))
        mkpath(dirname(target_path))
    end
    save(target_path, fig)
    if verbose
        @info "saved plot to " target_path
    end
end

function savePlot(
    fig, target_dir::String, target_fn::String; overwrite::Bool=false, verbose::Bool=true
)
    savePlot(fig, joinpath(target_dir, target_fn); overwrite, verbose)
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
    if lat < 0
        val = "$(abs(lat))°S"
    elseif lat == 0
        val = "0°"
    else
        val = "$(lat)°N"
    end
    return val
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
                xlabel_rotate = xlabel_rotate
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
                xlabel_rotate = xlabel_rotate
            )
        end
    end
    return fig
end


function gradColors(values::AbstractArray; name::Symbol=:thermal, rev=true)
    cmap = cgrad(name, rev=rev, 0:0.1:1)
    vmin, vmax = minimum(values), maximum(values)
    # normalization function
    norm(v) = (v - vmin) / (vmax - vmin)
    colors = [cmap[norm(v)] for v in values]
    return colors
end


function addMinorGrid!(ax, data_x::AbstractArray, data_y::AbstractArray; by = 0.5)
    ax.xminorticks = floor(minimum(data_x)) : by : ceil(maximum(data_x))
    ax.yminorticks = floor(minimum(data_y)) : by : ceil(maximum(data_y))
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true
    return nothing
end


"""
    splitColormapAtZero(range_min, range_max; n = 128)

Build a colormap where blue colors are used for values in [range_min, 0)
and red colors are used for values in [0, range_max], with the boundary
placed at the correct relative position (handles asymmetric ranges too).
"""
function splitColormapAtZero(range_min, range_max; n::Int = 128)
    if range_min >= 0 || range_max <= 0
        @warn "Note: range does not straddle 0."
        #return get(ColorSchemes.colorschemes[:berlin], range(0, 1; length = n))
        return reverse(get(ColorSchemes.colorschemes[:redblue], range(0, 1; length = n)))
    end
    below_colors = reverse(get(ColorSchemes.colorschemes[:Blues], range(0, 1; length = n)))
    above_colors = get(ColorSchemes.colorschemes[:Reds], range(0, 1; length = n))
    combined_colors = vcat(below_colors, above_colors)
    # stops must be strictly increasing; nudge the boundary slightly so both
    # sides get their own stop right at ratio_below_zero
    ratio_below_zero = (0 - range_min) / (range_max - range_min)
    stops = vcat(
        range(0, ratio_below_zero; length = n),
        range(ratio_below_zero, 1; length = n) .+ eps() .* (1:n)
    )
    stops = stops ./ stops[end]  # renormalize to [0, 1]
    return Makie.cgrad(combined_colors, stops)
end


function addColorBar(
    fig, hm;
    pos_legend::Union{Nothing, NamedTuple} = nothing,
    orient_legend::Symbol = :vertical,
    legend_label::String = "",
    fontsize::Int = 20,
    colorbar_size::Int = 15
)
    if orient_legend == :vertical
        cbgrid = GridLayout(3, 1)
        if pos_legend.x == 0
            fig[:, pos_legend.y] = cbgrid
        else
            fig[pos_legend.x, pos_legend.y] = cbgrid
        end
        rowsize!(cbgrid, 1, Relative(0.1))
        rowsize!(cbgrid, 2, Relative(0.8))
        rowsize!(cbgrid, 3, Relative(0.1))
        Colorbar(
            cbgrid[2, 1], hm, 
            width = Fixed(colorbar_size),
            labelsize = fontsize,
            ticklabelsize = fontsize - 2, 
            label = legend_label, 
            ticksvisible = false
        )
    else
        cbgrid = GridLayout(1, 3)
        if pos_legend.y == 0
            fig[pos_legend.x, :] = cbgrid
        else
            fig[pos_legend.x, pos_legend.y] = cbgrid
        end
        colsize!(cbgrid, 1, Relative(0.1))
        colsize!(cbgrid, 2, Relative(0.8))
        colsize!(cbgrid, 3, Relative(0.1))
        Colorbar(
            cbgrid[1,2], hm; 
            height = Fixed(colorbar_size),
            flipaxis = false,
            vertical = false,
            ticksvisible = false,
            labelsize = fontsize,
            ticklabelsize = fontsize - 2, 
            label = legend_label
        )
    end
end