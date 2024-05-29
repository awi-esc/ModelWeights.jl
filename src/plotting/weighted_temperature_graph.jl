using NetCDF
using Statistics
using CairoMakie

#CairoMakie.activate!()

include(joinpath("..", "load-data-utils.jl"));
include("plot-utils.jl")

workDir = joinpath(PATH_TO_WORK_DIR, "weighted_temperature_graph", "weighted_temperature_graph");
targetDir = joinpath("plots-replicated-with-julia", "weighted_temperature_graph");
mkpath(targetDir)

path_to_data = joinpath(workDir, "temperature_anomalies.nc")

# use ncdump to see the variables contained within the .nc-file.
#  or show contents like so:
#ncfile = NetCDF.open(path_to_data)
temp_anomalies = ncread(path_to_data, "tas");
size(temp_anomalies)

time =  ncread(path_to_data, "time");
size(time)

model_ensemble =  ncread(path_to_data, "model_ensemble");
size(model_ensemble)


unweighted_avg = mean(temp_anomalies, dims = 2);
size(unweighted_avg)

xticks = 1960:20:2100
timeInYears = trunc.(Int, 1960 .+ (time ./ 365))

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
    range66 = [quantile(row, [0.167, 0.833]) for row in eachrow(temp_anomalies)];
    lower =  unweighted_avg .- getindex.(range66, 1)
    upper = unweighted_avg .+ getindex.(range66, 2)
    band!(ax, timeInYears, vec(lower), vec(upper), color = (:red, 0.2))
    
    points = Point2f.(timeInYears, vec(unweighted_avg))
    
    for col in 1:size(temp_anomalies)[2]
        y = temp_anomalies[:, col]  # Extract the column
        lines!(ax, timeInYears, y, color = :gray70, label = "ensemble members")
    end
    lines!(ax, timeInYears, vec(unweighted_avg), color = :red, label = "Non-weighted mean")
    
    axislegend(ax, merge = true, position = :lt)
end
f

save(joinpath(targetDir, getCurrentTime() * "_weighted_temperature_graph.png"), f);

