""" plotValsOnMap!(fig::Figure, means::AbstractArray, title::String;
                    colors=nothing, color_range=nothing, high_clip=(1,0,0),
                    low_clip=(0,0,1), pos=(x=1, y=1), pos_legend=nothing
                    )
    
Plot contours of world with an overlayed heatmap of the input data.

# Arguments:
- `pos::NamedTuple(x::Int,y::Int)`: position of plot in `fig`
- `color_range::Union{Nothing, Tuple}`: Real values, outside this range, a single color is used.
"""
function plotValsOnMap!(
    fig::Figure,
    means::AbstractArray,
    title::String;
    colors = nothing,
    color_range::Union{Nothing, Tuple} = nothing,
    pos::NamedTuple = (x = 1, y = 1),
    pos_legend::Union{Nothing, NamedTuple} = (x = 1, y = 2),
    orient_legend::Symbol = :vertical,
    xlabel::String = "Longitude",
    ylabel::String = "Latitude",
    xlabel_rotate::Number = 0,
    nb_ticks::Union{Int,Nothing} = nothing,
    east_west_labels::Bool = false,
    alpha::Number = 0.8
)
    means = Data.sortLongitudesWest2East(means)
    dims_lat = Array(dims(means, :lat))
    dims_lon = Array(dims(means, :lon))
    dims_lon = Data.lon360to180.(dims_lon)

    # scaling plot 
    lon_min, lon_max = minimum(dims_lon) - 1, maximum(dims_lon) + 1
    lat_min, lat_max = minimum(dims_lat) - 1, maximum(dims_lat) + 1
    lon = range(lon_min, stop = lon_max, length = length(dims_lon))
    lat = range(lat_min, stop = lat_max, length = length(dims_lat))

    # axis ticks and labels
    xticks = isnothing(nb_ticks) ? [ceil(dims_lon[1]), 0, round(dims_lon[end])] : dims_lon
    yticks = isnothing(nb_ticks) ? [ceil(dims_lat[1]), 0, round(dims_lat[end])] : dims_lat
    lon_labels = east_west_labels ? longitude2EastWest.(xticks) : string.(xticks)
    lat_labels = east_west_labels ? latitude2NorthSouth.(yticks) : string.(yticks)

    step_lon = isnothing(nb_ticks) ? 1 : Int(round(length(lon_labels) / nb_ticks))
    step_lat = isnothing(nb_ticks) ? 1 : Int(round(length(lat_labels) / nb_ticks))
    x_ticks_labels = (xticks[1:step_lon:end], lon_labels[1:step_lon:end])
    y_ticks_labels = (yticks[1:step_lat:end], lat_labels[1:step_lat:end])

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
    hm = isnothing(color_range) ?
        heatmap!(lon, lat, Array(means); colormap=colors, alpha=alpha) :
        heatmap!(
            lon,
            lat,
            Array(means);
            colormap = colors,
            alpha = alpha,
            colorrange = color_range,
            highclip = colors[end],
            lowclip = colors[1]
        )
    lines!(GeoMakie.coastlines(); color = :black)
    if !isnothing(pos_legend)
        if orient_legend == :vertical
            Colorbar(fig[pos_legend.x, pos_legend.y], hm, height = Relative(2 / 3))
        else
            Colorbar(fig[pos_legend.x, pos_legend.y], hm, width = Relative(2 / 3), vertical = false)
        end
        # the width argument is relative to the width of the column#, width=Relative(0.1));
    end
    return nothing
end

function plotValsOnMap(    
    means::AbstractArray,
    title::String;
    colors = nothing,
    color_range::Union{Nothing, Tuple} = nothing,
    pos::NamedTuple = (x = 1, y = 1),
    pos_legend::Union{Nothing, NamedTuple} = (x = 1, y = 2),
    orient_legend::Symbol = :vertical,
    xlabel::String = "Longitude",
    ylabel::String = "Latitude",
    xlabel_rotate::Number = 0,
    nb_ticks::Union{Int,Nothing} = nothing,
    east_west_labels::Bool = false,
    alpha::Number = 0.8
)
    f = Figure()
    plotValsOnMap!(
        f, means, title; 
        colors, color_range, 
        pos, pos_legend, orient_legend, 
        xlabel, ylabel, xlabel_rotate, nb_ticks, east_west_labels,
        alpha
    )
    return f
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

    unit = get(data.properties, "units", "")
    xlab = isempty(unit) ? "" : "Transport in " * unit
    ax = Axis(
        fig[1, 1],
        title = join(["AMOC strength (at 26.5°N)", t], "\n"),
        xlabel = xlab
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
    plotEnsembleSpread(
        data::YAXArray; fn::Function=Statistics.median, title::String="", ylabel::String=""
    )

Plot models against values of all model members. Add summary statistics (default: median) 
when computed across all model members and when average values across members of each model 
are computed first.

# Arguments:
- `data::YAXArray`: must have single dimension :member.
- `fn::Function=Statistics.median`: function to summarize data.
- `title::String`: (optional) plot title.
- `ylabel::String`: (optional) label for y-axis.
"""
function plotEnsembleSpread(
    data::YAXArray; fn::Function=Statistics.median, title::String="", ylabel::String=""
)
    if !hasdim(data, :member) || length(otherdims(data, :member)) >= 1
        msg = "To plot ensemble spread input data must have a single dimension :member"
        throw(ArgumentError(msg))
    end
    members = collect(lookup(data, :member))
    df = Data.setLookupsFromMemberToModel(data, ["member"]) 
    models = collect(lookup(df, :model))
    model_labels = unique(models)
    model_ints  = collect(1:length(members))
    for (i, m) in enumerate(model_labels)
        indices = findall(models .== m)
        model_ints[indices] .= i
    end

    counts = Data.countMap(model_ints)
    model_labels = map(((i, x),) -> x * " ($(string(counts[i])))", enumerate(model_labels))

    units = get(df.properties, "units", nothing)
    if !isnothing(units)
        if length(unique(units)) != 1
            @warn "Not all data is defined in the same units! Found: $(units)!"
        end
    end

    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "Models",
        ylabel = ylabel,
        title = title,
        xticks = (1:length(model_labels), model_labels),
        xticklabelrotation = pi / 2
    )
    scatter!(ax, model_ints, collect(df.data))

    summarized_members = fn(df)
    summarized_avg_models = fn(Data.summarizeMembers(data))
    fn_name = string(nameof(fn)) 
    l1, l2 =  fn_name * " across all members", fn_name * " members first summarized within model"
    hlines!(ax, summarized_members; linestyle = :dash, color=:magenta, label=l1)
    hlines!(ax, summarized_avg_models; linestyle = :dash, color=:green, label=l2)
    axislegend(fontsize = 10)
    return fig
end

"""
    plotTimeseries!(ax::Axis, vals::AbstractArray;)

Plot timeseries of data vector `vals`.
"""
function plotTimeseries!(
    ax::Axis,
    vals::AbstractArray;
    uncertainties::Union{AbstractArray, Nothing} = nothing,
    color_line::Symbol = :darkred,
    color_unc::Symbol = :darkred,
    label::String = "",
    label_unc::String = "",
    linestyle::Symbol = :solid,
    linewidth = 3,
    alpha = 0.5
)
    plots = []
    timesteps = Array(dims(vals, :time))
    if typeof(timesteps[1]) == DateTime
        timesteps = map(x -> Dates.year(x), timesteps)
    end

    lineplot = lines!(
        ax,
        timesteps,
        vec(coalesce.(vals, NaN)),
        color = color_line, # color = (color_line, alpha), 
        label = label,
        linestyle = linestyle,
        linewidth = linewidth,
    )
    push!(plots, lineplot)
    if !isnothing(uncertainties)
        bandplot = band!(
            ax,
            timesteps,
            vec(coalesce.(uncertainties[confidence=At("lower")], NaN)),
            vec(coalesce.(uncertainties[confidence=At("upper")], NaN)),
            color = (color_unc, alpha),
            label = label_unc
        )
        push!(plots, bandplot)
    end
    return plots
end


"""
    plotTimeseries(ax::Axis, vals::AbstractArray;)

Plot timeseries of data vector `data`.

# Arguments:
- `data::YAXArray`: must have dimension 'time' and possibly one other dimension
"""
function plotTimeseries(
    data::YAXArray;
    linestyle::Symbol = :solid,
    linewidth = 3,
    xlabel = "time",
    ylabel = "",
    title = "",
    colors::AbstractArray=[],
    n_step::Int = 10
)
    timesteps = Array(dims(data, :time))
    nt = length(timesteps)
    timesteps = 1:nt

    if nt < n_step
        xticks = 1:nt
    else
        xticks = [1, (n_step : n_step : nt)...]
    end
    xtick_labels = string.(xticks)

    f = Figure(); 
    ax = Axis(f[1,1], xlabel=xlabel, ylabel=ylabel, title=title, xticks=(xticks, xtick_labels));

    nb_dims = ndims(data)
    if nb_dims > 2
        throw(ArgumentError("Timeseries can only be plotted for data with :time dimension and just one other dimension."))
    end
    if nb_dims == 1
        lines!(
            ax,
            timesteps,
            vec(coalesce.(data, NaN)),
            linestyle = linestyle,
            linewidth = linewidth
        )
    else
        idx_time_dim = Data.indexDim(data, :time)
        idx_other_dim = idx_time_dim == 1 ? 2 : 1
        n = size(data, idx_other_dim)
        plots = Vector(undef, n)
        for idx in eachindex(1:n)
            indices = idx_time_dim == 1 ? [:, idx] : [idx , :]
            if !isempty(colors)
                plots[idx] = lines!(
                    ax,
                    timesteps,
                    vec(coalesce.(data[indices...], NaN)),
                    linestyle = linestyle,
                    linewidth = linewidth,
                    color = colors[idx]
                )
            else
                plots[idx] = lines!(
                    ax,
                    timesteps,
                    vec(coalesce.(data[indices...], NaN)),
                    linestyle = linestyle,
                    linewidth = linewidth
                )
            end
        end
        Legend(
            f[2,1], plots, Array(dims(data)[idx_other_dim]), framevisible=false, 
            orientation=:horizontal, nbanks = div(n,4) +1, 
            labelsize=10
        )
    end
    return f
end


"""
    plotTempGraph

    # TODO: change quantileLabels and uncertainty_ranges, the latter should contain 
    # the quantile labels
"""
function plotTempGraph(
    data::YAXArray,
    averages::NamedTuple,
    uncertainty_ranges::NamedTuple,
    title::String;
    ylabel::String = "",
)
    f = Figure()
    years = Dates.year.(Array(dims(data, :time)))
    xticks = years[1]:20:years[end]
    units = unique(get(data.properties, "units", []))
    if length(units) > 1
        throw(ArgumentError("multiple units in data!"))
    elseif isempty(units)
        @warn "plotTempGraph: no units given in metadata!"
    end
    unit = !isempty(units) ? units[1] : ""
    ylab = isempty(ylabel) ? (isempty(unit) ? "Temperature" : "Temperature in $unit") : ylabel
    ax = Axis(
        f[1, 1],
        xticks = (xticks, string.(xticks)),
        limits = ((years[1] - 10, years[end] + 10), (-1, 5)),
        title = title,
        xlabel = "Year",
        ylabel = ylab,
    )
    # add ranges TODO: make label with fn argument
    lower_unw = collect(uncertainty_ranges.unweighted[confidence=At("lower")])
    upper_unw = collect(uncertainty_ranges.unweighted[confidence=At("upper")])
    band!(
        ax,
        years,
        coalesce.(lower_unw, NaN),
        coalesce.(upper_unw, NaN),
        color = (:red, 0.2),
        label = "Non-weighted 16.7-83.3 perc range",
    )
    lower_weighted = collect(uncertainty_ranges.weighted[confidence=At("lower")])
    upper_weighted = collect(uncertainty_ranges.weighted[confidence=At("upper")])
    band!(
        ax,
        years,
        coalesce.(lower_weighted, NaN),
        coalesce.(upper_weighted, NaN),
        color = (:green, 0.2),
        label = "Weighted 16.7-83.3perc range",
    )
    # add results for each model model seperately
    for m in dims(data, :member)
        y = data[member=At(m)]
        lines!(ax, years, Array(y), color = :gray80, linewidth=1, label = "Ensemble members")
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


function plotECDF(vals::Vector{<:Number}; variable::String="", title::String="Empirical CDF")
    fn = ecdf(vals)
    indices = sortperm(vals)
    sorted = unique(vals[indices]) 
    f = Figure()
    ax = Axis(f[1,1], xticks=vals, title = title, ylabel = "Cumulative Probability", xlabel=variable)
    ylims!(0,1)
    for (i, v) in enumerate(sorted[1:end-1])
        p = fn(v)
        next_v = sorted[i + 1]
        scatter!(v, p, color=:grey)
        lines!([v, sorted[i+1]], [p, p], color=:grey)
        lines!([next_v, next_v], [p, fn(next_v)], color=:grey)
    end
    return f
end


function plotPDF(xs, ys, xlabel::String;
    samples::Vector{<:Number},
    selected_samples::Vector{<:Number},
    selected_names::Vector{String},
    label::String="",
    label_samples::String="",
    markersize::Number=20,
    ncols_leg::Number=2,
    label_size::Number=12,
    colors
)
    f_pdf = Figure()
    ax = Axis(f_pdf[1,1], xlabel = xlabel, ylabel = "Density")
    if isempty(label)
        scatter!(ax, xs, ys)
    else
        scatter!(ax, xs, ys, label=label)
    end

    # also plot samples if given
    scatter_y = rand(Distributions.Normal(0, 0.02), length(samples))
    if isempty(label_samples)
        scatter!(ax, samples, scatter_y)
    else
        scatter!(ax, samples, scatter_y, label=label_samples)
    end
    # label and color models where the ecs based performance weight deviates much from the historical based
    for (i, s) in enumerate(selected_samples)   
        scatter!(ax, s, scatter_y[i], label="$(selected_names[i])", color=colors[i], markersize=markersize)
    end
    axislegend(labelsize=label_size, nbanks=ncols_leg)
    return f_pdf
end