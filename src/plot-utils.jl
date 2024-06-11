using ColorSchemes
using Dates

using CairoMakie    
CairoMakie.activate!(type = "svg")


function plotDistMatrices(distMat, climateVar, models, modelRefs)
    size_inches = (8, 5)
    size_pt = 72 .* size_inches
    fig=Figure(size= size_pt, fontsize=12)

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


function plotPerformanceMetric(data, climateVar, models)
    size_inches = (6.5, 6)
    size_pt = 72 .* size_inches
    fig=Figure(size= size_pt, fontsize=12)

    xs = 1:length(models)
    ax = Axis(fig[1,1], 
        xticks = (xs, models), 
        xticklabelrotation = pi/4,
        title = "RMS error for " * climateVar,
        xlabel = "Model", 
        ylabel = "RMS error " * climateVar
    );
    barplot!(ax, xs, data)
    return fig
end

function longitude2EastWest(lon)
    return lon > 0 ? "$(lon)°E" : "$(abs(lon))°W"
end

function latitude2NorthSouth(lat)
    return lat < 0 ? "$(lat)°S" : "$(abs(lat))°N"
end


function getCurrentTime()
    currentDay = string(today()) * '_';
    currentTime = Dates.format(now(), "HH_MM");
    return currentDay * currentTime
end