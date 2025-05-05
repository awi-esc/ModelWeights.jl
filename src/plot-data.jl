using ColorSchemes
using Dates

""" plotValsOnMap!(fig::Figure, means::AbstractArray, title::String;
                    colors=nothing, color_range=nothing, high_clip=(1,0,0),
                    low_clip=(0,0,1), pos=(x=1, y=1), pos_legend=nothing
                    )
    
Plot contours of world with an overlayed heatmap of the input data.

# Arguments:
- `nb_ticks`: if nothing (default), just min lat/lon, max lat/lon and 0 are labeled.
"""
function plotValsOnMap!(
    fig::Figure,
    means::AbstractArray,
    title::String;
    colors = nothing,
    color_range = nothing,
    high_clip = (1, 0, 0),
    low_clip = (0, 0, 1),
    pos = (x = 1, y = 1),
    pos_legend = (x = 1, y = 2),
    xlabel = "Longitude",
    ylabel = "Latitude",
    xlabel_rotate = pi / 4,
    nb_ticks::Union{Int,Nothing} = nothing,
    east_west_labels = false,
)
    means = sortLongitudesWest2East(means)
    dims_lat = Array(dims(means, :lat))
    dims_lon = Array(dims(means, :lon))
    if any(x -> x > 179, dims_lon)
        dims_lon = lon360to180.(dims_lon)
    end

    # scaling plot 
    lon_min, lon_max = minimum(dims_lon) - 1, maximum(dims_lon) + 1
    lat_min, lat_max = minimum(dims_lat) - 1, maximum(dims_lat) + 1
    lon = range(lon_min, stop = lon_max, length = length(dims_lon))
    lat = range(lat_min, stop = lat_max, length = length(dims_lat))

    # axis ticks and labels
    xticks = isnothing(nb_ticks) ? [ceil(dims_lon[1]), 0, round(dims_lon[end])] : dims_lon
    yticks = isnothing(nb_ticks) ? [ceil(dims_lat[1]), 0, round(dims_lat[end])] : dims_lat
    lonLabels = east_west_labels ? longitude2EastWest.(xticks) : string.(xticks)
    latLabels = east_west_labels ? latitude2NorthSouth.(yticks) : string.(yticks)

    step_lon = isnothing(nb_ticks) ? 1 : Int(round(length(lonLabels) / nb_ticks))
    step_lat = isnothing(nb_ticks) ? 1 : Int(round(length(latLabels) / nb_ticks))
    x_ticks_labels = (xticks[1:step_lon:end], lonLabels[1:step_lon:end])
    y_ticks_labels = (yticks[1:step_lat:end], latLabels[1:step_lat:end])

    ax = Axis(
        fig[pos.x, pos.y],
        title = title,
        xlabel = xlabel,
        ylabel = ylabel,
        xticklabelrotation = xlabel_rotate,
        xticks = x_ticks_labels,
        yticks = y_ticks_labels,
        limits = ((lon_min, lon_max), (lat_min, lat_max)),
    )
    if isnothing(colors)
        colors = reverse(ColorSchemes.redblue.colors)
    end
    hm =
        isnothing(color_range) ?
        heatmap!(ax, lon, lat, Array(means), colormap = colors, alpha = 0.8) :
        heatmap!(
            ax,
            lon,
            lat,
            Array(means),
            colormap = colors,
            alpha = 0.8,
            colorrange = color_range,
            highclip = high_clip,
            lowclip = low_clip,
        )
    lines!(GeoMakie.coastlines(); color = :black)
    if !isnothing(pos_legend)
        Colorbar(fig[pos_legend.x, pos_legend.y], hm, height = Relative(2 / 3))
        # the width argument is relative to the width of the column#, width=Relative(0.1));
    end
    return nothing
end


""" plotHistAtPos(data::AbstractArray, location::Dict)

Plot histogram of all data for a specific `location`.

# Arguments:
- `data`: dimensions 'lon' (from -180° to 180°), 'lat' (-90° to 90°)
- `location`: must have keys 'name', 'lon', 'lat'
"""
function plotHistAtPos(data::YAXArray, location::Dict)
    longitudes = Array(dims(data, :lon))
    longitudes = ifelse(any(longitudes .> 180), lon360to180.(longitudes), longitudes)
    if location["lon"] > 180
        location["lon"] = lon360to180(location["lon"])
    end
    coords = getClosestGridPoint(location, longitudes, Array(dims(data, :lat)))
    # no idea why At doesnt work for negative values at longitude dimension...
    # At would drop dimensions directly..
    data_loc =
        data[lat=Where(elem->elem==coords["lat"]), lon=Where(elem->elem==coords["lon"])]
    data_loc = dropdims(data_loc, dims = :lon)
    data_loc = dropdims(data_loc, dims = :lat)

    grid_lat = latitude2NorthSouth(coords["lat"])
    grid_lon = longitude2EastWest(coords["lon"])
    t1 =
        "Variable: " *
        data.properties["variable_id"] *
        " Experiment: " *
        data.properties["experiment_id"]
    t2 = "near " * location["name"] * "(" * grid_lat * "," * grid_lon * ")"

    fig = getFigure((14, 10), 16)
    ax = Axis(fig[1, 1], title = join([t1, t2], "\n"), xlabel = data.properties["units"])
    hist!(ax, vec(data_loc))
    return fig
end


""" plotAMOC(data::DimArray)

Plot AMOC strength for variable "amoc".
"""
function plotAMOC(data::YAXArray)
    fig = getFigure((8, 5), 12)
    t = "variable: amoc, experiment: " * data.properties["experiment_id"]

    unit = data.properties["units"]
    ax = Axis(
        fig[1, 1],
        title = join(["AMOC strength (at 26.5°N)", t], "\n"),
        xlabel = "Transport in " * unit,
    )
    if length(dims(data)) == 1
        # data for full period
        hist!(ax, Array(data))
    else
        # seasonal data
        for i = 1:4
            hist!(ax, Array(data[i, :]), scale_to = -0.6, offset = i, direction = :x)
        end
        ax.xlabel = "season"
    end
    return fig
end


"""
    plotEnsembleSpread(data::AbstractArray, lon::Number, lat::Number)

Create figure with boxplots for each model in `data` that have several ensemble members.
"""
function plotEnsembleSpread(data::YAXArray, lon::Number, lat::Number)
    # models = unique(dims(data, :model));
    # models_ensembles = filter(x -> length(dims(data[model = Where(m -> m == x)], :model)) > 1, models);
    # data_ensembles = data[model = Where(x -> x in models_ensembles)];
    data_ensembles = data
    models_ensembles = unique(dims(data, :model))
    # translate list of unique models into list of integers for boxplot
    categories = Array(dims(data_ensembles, :model))
    categoriesInts = collect(1:length(categories))
    for (i, m) in enumerate(models_ensembles)
        categoriesInts[findall(categories .== m)] .= i
    end

    fig = getFigure((16, 8), 18)
    t1 = data.properties["long_name"] * " (" * data.properties["variable_id"] * ")"
    t = join(["at ", longitude2EastWest(lon), latitude2NorthSouth(lat)], " ")
    t2 = "Spread of models with several ensemble members " * t

    ax = Axis(
        fig[1, 1],
        xlabel = "Models",
        ylabel = data.properties["units"],
        title = join([t1, t2], "\n"),
        xticks = (collect(1:length(models_ensembles)), models_ensembles),
        xticklabelrotation = pi / 2,
    )

    if !(lon in Array(dims(data_ensembles, :lon)))
        msg = "location and data are not given in the same longitudes scale, \
         use either (0° to 360°) or (-180° to 180°)!"
        throw(ArgumentError(msg))
    end
    values =
        dropdims(data_ensembles[lon=Where(x->x==lon), lat=Where(x->x==lat)], dims = :lon)
    boxplot!(ax, categoriesInts, Array(dropdims(values, dims = :lat)))
    return fig
end


function plotTimeseries!(
    ax::Axis,
    vals::AbstractArray;
    uncertainties::Union{AbstractArray,Nothing} = nothing,
    color_line::Symbol = :darkred,
    color_range::Symbol = :red,
    label::String = "",
    label_unc::String = "",
    linestyle::Symbol = :solid,
    linewidth = 2,
    alpha = 0.2,
)
    timesteps = Array(dims(vals, :time))
    if typeof(timesteps[1]) == DateTime
        timesteps = map(x -> Dates.year(x), timesteps)
    end

    lines!(
        ax,
        timesteps,
        vec(coalesce.(vals, NaN)),
        color = color_line, # color = (color_line, alpha), 
        label = label,
        linestyle = linestyle,
        linewidth = linewidth,
    )
    if !isnothing(uncertainties)
        band!(
            ax,
            timesteps,
            vec(coalesce.(uncertainties[confidence=At("lower")], NaN)),
            vec(coalesce.(uncertainties[confidence=At("upper")], NaN)),
            color = (color_range, alpha),
            label = label_unc,
        )
    end
    return nothing
end


"""
    plotTempGraph

    # TODO: change quantileLabels and uncertaintyRanges, the latter should contain 
    # the quantile labels
"""
function plotTempGraph(
    data::YAXArray,
    averages::NamedTuple,
    uncertaintyRanges::NamedTuple,
    title::String;
    ylabel::String = "",
)
    f = Figure()
    years = Dates.year.(Array(dims(data, :time)))
    xticks = years[1]:20:years[end]
    ax = Axis(
        f[1, 1],
        xticks = (xticks, string.(xticks)),
        limits = ((years[1] - 10, years[end] + 10), (-1, 5)),
        title = title,#"Temperature anomaly relative to " * name_ref_period,
        xlabel = "Year",
        ylabel = isempty(ylabel) ? "Temperature in " * data.properties["units"] : ylabel,
    )
    # add ranges TODO: make label with fn argument
    lowerUnw = uncertaintyRanges.unweighted[confidence=At("lower")]
    upperUnw = uncertaintyRanges.unweighted[confidence=At("upper")]
    band!(
        ax,
        years,
        vec(lowerUnw),
        vec(upperUnw),
        color = (:red, 0.2),
        label = "Non-weighted 16.7-83.3 perc range",
    )

    lowerWeighted = uncertaintyRanges.weighted[confidence=At("lower")]
    upperWeighted = uncertaintyRanges.weighted[confidence=At("upper")]
    band!(
        ax,
        years,
        vec(lowerWeighted),
        vec(upperWeighted),
        color = (:green, 0.2),
        label = "Weighted 16.7-83.3perc range",
    )

    # add results for each model model seperately
    for m in dims(data, :member)
        y = data[member=At(m)]
        lines!(ax, years, Array(y), color = :gray80, label = "Ensemble members")
    end

    lines!(ax, years, vec(averages.unweighted), color = :red, label = "Non-weighted mean")
    lines!(ax, years, vec(averages.weighted), color = :green, label = "Weighted mean")
    axislegend(ax, merge = true, position = :lt)
    return f
end


function makeScatterPlot(
    xs::AbstractArray{<:Union{Missing,Number}},
    ys::AbstractArray{<:Union{Missing,Number}};
    captions::NamedTuple = (x = "", y = "", title = ""),
    xtick_labels::Union{Vector{String},Nothing} = nothing,
    xticklabelrotation::Number = pi / 2,
    legend::NamedTuple = (label = "", position = :rc, color = :red),
    greyed_area::NamedTuple = (
        y1 = Inf,
        y2 = Inf,
        label = "",
        summary_val = Inf,
        summary_stat = "",
    ),
    add_lines::Bool = true,
)
    f = Figure()
    ax =
        isnothing(xtick_labels) ?
        Axis(f[1, 1], title = captions.title, xlabel = captions.x, ylabel = captions.y) :
        Axis(
            f[1, 1],
            xticks = (xs, xtick_labels),
            title = captions.title,
            xlabel = captions.x,
            ylabel = captions.y,
            xticklabelrotation = xticklabelrotation,
        )
    if add_lines
        lines!(ax, Array(xs), Array(ys), color = legend.color, label = legend.label)
    end
    scatter!(ax, Array(xs), Array(ys), color = legend.color, label = legend.label)

    if !isinf(greyed_area.y1) && !isinf(greyed_area.y2)
        band!(
            xs,
            greyed_area.y1,
            greyed_area.y2,
            color = (:gray, 0.5),
            label = greyed_area.label,
        )
        lines!(
            ax,
            xs,
            repeat([greyed_area.summary_val], length(xs)),
            color = (:gray, 0.5),
            label = greyed_area.summary_stat,
        )
    end
    if !isempty(legend.label)
        axislegend(ax, merge = true, position = legend.position)
    end
    return (f, ax)
end
