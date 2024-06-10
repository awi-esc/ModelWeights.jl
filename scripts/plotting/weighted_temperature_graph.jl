using NetCDF
using Statistics
using CairoMakie
using Interpolations

include(joinpath(@__DIR__, "..", "..", "src", "load-data-utils.jl"));
include(joinpath(@__DIR__, "..", "..", "src", "plot-utils.jl"));

workDir = joinpath(PATH_TO_WORK_DIR, "weighted_temperature_graph", "weighted_temperature_graph");
targetDir = joinpath("plots-replicated-with-julia", "weighted_temperature_graph");
mkpath(targetDir)

path2Data = joinpath(workDir, "temperature_anomalies.nc")

# use ncdump to see the variables contained within the .nc-file.
#  or show contents like so:
#ncfile = NetCDF.open(path_to_data)
tempAnomalies = ncread(path2Data, "tas");
size(tempAnomalies)

time =  ncread(path2Data, "time");
size(time)
modelEnsemble =  ncread(path2Data, "model_ensemble");
size(modelEnsemble)
nbModels = length(modelEnsemble);

xticks = 1960:20:2100
timeInYears = trunc.(Int, 1960 .+ (time ./ 365));

function getInterpolatedWeightedQuantiles(quantiles, vals, weights=nothing)
    if isnothing(weights)
        weights = ones(length(vals));
    end
    indicesSorted = sortperm(vals);
    weightsSorted = weights[indicesSorted];
    weightedQuantiles = cumsum(weightsSorted) - 0.5 * weightsSorted;
    weightedQuantiles = reshape(weightedQuantiles, length(weightedQuantiles), 1);
    weightedQuantiles = (weightedQuantiles .- minimum(weightedQuantiles)) ./ maximum(weightedQuantiles);
    
    interp_linear = Interpolations.linear_interpolation(vec(weightedQuantiles), 
                                                        vals[indicesSorted],
                                                        extrapolation_bc=Line()
                                                        );
     
    return interp_linear(quantiles)
end

weights = ncread(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "weights.nc"),
                 "weight");

quantiles = [0.167, 0.833];
unweightedRanges = [];
weightedRanges = [];
for i_row in axes(tempAnomalies, 1)
    lower, upper = getInterpolatedWeightedQuantiles(quantiles, tempAnomalies[i_row,:], weights);
    push!(weightedRanges, [lower, upper]);
    lower, upper = getInterpolatedWeightedQuantiles(quantiles, tempAnomalies[i_row,:]);
    push!(unweightedRanges, [lower, upper]);
end

begin 
    f = Figure(); 
    ax = Axis(f[1,1], 
        xticks = (xticks, string.(xticks)), 
        limits = ((minimum(timeInYears)-10, maximum(timeInYears)+10), (-1, 5)),
        title = "Temperature anomaly relative to 1981-2010",
        xlabel = "Year", 
        ylabel = "Temperature anomaly °C"
    );
    # add ranges 
    lowerUnw = getindex.(unweightedRanges, 1);
    upperUnw = getindex.(unweightedRanges, 2);
    band!(ax, timeInYears, vec(lowerUnw), vec(upperUnw), color = (:red, 0.2), 
          label = "Non-weighted 16.7-83.3perc range");
    
    lowerWeighted = getindex.(weightedRanges, 1);
    upperWeighted = getindex.(weightedRanges, 2);
    band!(ax, timeInYears, vec(lowerWeighted), vec(upperWeighted), color = (:green, 0.2), 
          label = "Weighted 16.7-83.3perc range");
        
    # add results for each model model seperately
    for i_model in 1:nbModels
        y = tempAnomalies[:, i_model] 
        lines!(ax, timeInYears, y, color = :gray80, label = "ensemble members")
    end
    unweightedAvg =  mean(tempAnomalies, dims = 2);
    weightedAvg = sum(repeat(weights', size(tempAnomalies,1), 1) .* tempAnomalies, dims=2);
    lines!(ax, timeInYears, vec(unweightedAvg), color = :red, label = "Non-weighted mean")
    lines!(ax, timeInYears, vec(weightedAvg), color = :green, label = "Weighted mean")


    axislegend(ax, merge = true, position = :lt)
end
f

save(joinpath(targetDir, getCurrentTime() * "_weighted_temperature_graph.png"), f);

