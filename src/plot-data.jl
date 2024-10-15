using ColorSchemes

""" plotMeansOnMap(means::DimArray, title::String, colors=nothing)
    
Plot contours of world with an overlayed heatmap that shows the data which 
correspond to mean value for each position in considered grid. 
"""
function plotMeansOnMap(means::DimArray, title::String, colors=nothing)
    means = sortLongitudesWest2East(means);
    dims_lat = Array(dims(means, :lat));
    dims_lon = Array(dims(means, :lon));
    if any(x -> x > 179, dims_lon)
        dims_lon = lon360to180.(dims_lon);
    end
    
    # scaling plot 
    lon_min, lon_max = minimum(dims_lon) -1, maximum(dims_lon) + 1;
    lat_min, lat_max = minimum(dims_lat) -1, maximum(dims_lat) + 1;
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
    if isnothing(colors)
        colors = reverse(ColorSchemes.redblue.colors)
    end
    hm = heatmap!(ax, lon, lat, Array(means), colormap=colors, alpha=0.8);
    lines!(GeoMakie.coastlines(); color=:black);
    Colorbar(fig[1,2], hm);

    return fig
end


""" plotHistAtPos(data::DimArray, location::Dict)

Plot histogram of all data for a specific `location`.

# Arguments:
- `data::DimArray`: dimensions 'lon' (from -180° to 180°), 'lat' (-90° to 90°)
- `location::Dict`: keys 'name', 'lon', 'lat'
"""
function plotHistAtPos(data::DimArray, location::Dict)
    longitudes = Array(dims(data, :lon));
    longitudes = ifelse(any(longitudes .> 180), lon360to180.(longitudes), longitudes);
    if location["lon"] > 180
        location["lon"] = lon360to180(location["lon"]);
    end
    coords = getClosestGridPoint(location, longitudes, Array(dims(data, :lat)));
    # no idea why At doesnt work for negative values at longitude dimension...
    # At would drop dimensions directly..
    data_loc = data[lat=Where(elem-> elem == coords["lat"]), 
                    lon=Where(elem -> elem == coords["lon"])];
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


"""
    plotEnsembleSpread(data::DimArray, lon::Number, lat::Number)

Create figure with boxplots for each model in `data` that have several ensemble members.
"""
function plotEnsembleSpread(data::DimArray, lon::Number, lat::Number)
    # models = unique(dims(data, :model));
    # models_ensembles = filter(x -> length(dims(data[model = Where(m -> m == x)], :model)) > 1, models);
    # data_ensembles = data[model = Where(x -> x in models_ensembles)];
    data_ensembles = data;
    models_ensembles = unique(dims(data, :model))
    # translate list of unique models into list of integers for boxplot
    categories = Array(dims(data_ensembles, :model));
    categoriesInts = collect(1 : length(categories));
    for (i, m) in enumerate(models_ensembles)
        categoriesInts[findall(categories .== m)] .= i;
    end
    
    fig =  getFigure((16,8), 18);
    t1 = data.metadata["long_name"] * " (" * data.metadata["variable_id"] * ")";    
    t = join(["at ", longitude2EastWest(lon), latitude2NorthSouth(lat)], " ");
    t2 = "Spread of models with several ensemble members " * t; 

    ax = Axis(fig[1,1], 
              xlabel = "Models", 
              ylabel = data.metadata["units"], 
              title = join([t1, t2], "\n"), 
              xticks = (collect(1:length(models_ensembles)), models_ensembles), 
              xticklabelrotation = pi/2);
    
    if !(lon in Array(dims(data_ensembles, :lon)))
        msg = "location and data are not given in the same longitudes scale, \
         use either (0° to 360°) or (-180° to 180°)!"
        throw(ArgumentError(msg));
    end
    values = dropdims(
        data_ensembles[lon = Where(x -> x == lon),
        lat = Where(x -> x == lat)], 
        dims =:lon
    );
    boxplot!(ax, categoriesInts, Array(dropdims(values, dims=:lat)));
    return fig
end


function plotTempGraph(
    data::DimArray, 
    averages::NamedTuple, 
    uncertaintyRanges::NamedTuple,
    title::String,
    quantilesLabel::String
)
    f = Figure(); 
    years = Dates.year.(Array(dims(data, :time)))
    xticks = years[1] : 20 : years[end]
    ax = Axis(f[1,1], 
        xticks = (xticks, string.(xticks)), 
        limits = ((years[1]-10, years[end]+10), (-1, 5)),
        title = title,#"Temperature anomaly relative to " * name_ref_period,
        xlabel = "Year", 
        ylabel = "Temperature in " * data.metadata["units"]# "Temperature anomaly °C"
    );
    # add ranges TODO: make label with fn argument
    lowerUnw = getindex.(uncertaintyRanges.unweighted, 1);
    upperUnw = getindex.(uncertaintyRanges.unweighted, 2);
    band!(ax, years, vec(lowerUnw), vec(upperUnw), color = (:red, 0.2), 
         label = "Non-weighted 16.7-83.3 perc range");
    
    lowerWeighted = getindex.(uncertaintyRanges.weighted, 1);
    upperWeighted = getindex.(uncertaintyRanges.weighted, 2);
    band!(ax, years, vec(lowerWeighted), vec(upperWeighted), color = (:green, 0.2), 
          label = "Weighted 16.7-83.3perc range");
        
    # add results for each model model seperately
    for m in dims(data, :model)
        y = data[model = At(m)] 
        lines!(ax, years, Array(y), color = :gray80, label = "ensemble members")
    end
    
    lines!(ax, years, vec(averages.unweighted), color = :red, label = "Non-weighted mean")
    lines!(ax, years, vec(averages.weighted), color = :green, label = "Weighted mean")
    axislegend(ax, merge = true, position = :lt)
    return f
end


