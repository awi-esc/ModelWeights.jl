using Colors: Colorant

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
    # high_clip = (1, 0, 0),
    # low_clip = (0, 0, 1),
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
        #pos_legend = idx_plot == nb_subplots ? (x = 1:row, y = grid.ncols + 1) : nothing
        model = models[idx_plot]
        if hasdim(data, :member)
            plotValsOnMap!(
                fig[pos.x, pos.y],
                data[member=At(model)],
                "$model";
                # colors = colors,
                # high_clip = high_clip,
                # low_clip = low_clip,
                color_range = color_range_limits,
                #pos = pos,
                #pos_legend = pos_legend,
                xlabel = xlabel,
                ylabel = ylabel,
                xlabel_rotate = xlabel_rotate
            )
        else
            plotValsOnMap!(
                fig[pos.x, pos.y],
                data[model=At(model)],
                "$model";
                #colors = colors,
                # high_clip = high_clip,
                # low_clip = low_clip,
                color_range = color_range_limits,
                #pos = pos,
                #pos_legend = pos_legend,
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


function mapValsToColorScheme(vals; rev=true, scale=nothing)
    vmin, vmax = extrema(vals)
    color = vmin >= 0 ? :Reds : (vmax <= 0 ? :Blues : :RdBu)
    return cgrad(color; rev=rev, scale=scale)
end


function mapValsToColors(values::AbstractArray; colormap::Union{Nothing, Symbol} = nothing, scale=nothing)
    vmin, vmax = extrema(values)
    len_range = vmax - vmin
    normalized = (values .- vmin) ./ len_range
    if isnothing(colormap)
        colormap = vmin >= 0 ? :Reds : (vmax <= 0 ? :Blues : :RdBu)
    end
    return get(cgrad(colormap; scale), normalized)
end


# function getColormap(
#     values::AbstractArray, split_at::Number;
#     col_name::Symbol = :RdBu, 
#     scale::Union{Symbol, Nothing} = nothing,
#     rev::Bool = true
# )
#     vmin, vmax = extrema(values)
#     len_range = vmax - vmin
#     if !(split_at > vmin && split_at < vmax)
#         error("Cant split colors at $(split_at) as it is not in range of given values!")
#     end
#     frac_upto_x = (split_at - vmin) / len_range
#     return cgrad(col_name, [frac_upto_x]; scale, rev) # at this point btw. 0 and 1, colors are positioned 
# end



function _colorbarPosition(gp::GridPosition, pos::Symbol; use_span::Bool = false)
    row = first(gp.span.rows)
    col = first(gp.span.cols)
    if pos == :r
        c = col + 1
        return use_span ? gp.layout[1:row, c] : gp.layout[row, c]
    elseif pos == :l
        c = col - 1
        return use_span ? gp.layout[1:row, c] : gp.layout[row, c]
    elseif pos == :t
        r = row - 1
        return use_span ? gp.layout[r, 1:col] : gp.layout[r, col] 
    elseif pos == :b
        r = row + 1
        return use_span ? gp.layout[r, 1:col] : gp.layout[r, col]
    else
        error("'pos' must be one of :r, :l, :t, :b")
    end
end


function addColorBar(
    gp::GridPosition, 
    color_map,
    color_range, 
    pos::Symbol;
    legend_label::String = "",
    fontsize::Int = 15,
    colorbar_size::Int = 15,
    clip_vals_colorbar::Bool = false
)
    clip_kwargs = clip_vals_colorbar ? (;highclip = color_map[end], lowclip = color_map[1]) : (;)
    orientation_kwargs = pos in [:r, :l] ? 
        (; width = Fixed(colorbar_size), flipaxis = pos == :r, vertical = true) :
        (; height = Fixed(colorbar_size), flipaxis = pos == :t, vertical = false)
    
    return Colorbar(
        gp;
        colormap = color_map,
        colorrange = color_range,
        labelsize = fontsize,
        ticklabelsize = fontsize - 2, 
        label = legend_label, 
        ticksvisible = false,
        clip_kwargs...,
        orientation_kwargs...
    )
end