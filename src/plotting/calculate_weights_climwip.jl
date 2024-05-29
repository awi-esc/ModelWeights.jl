using CairoMakie

include(joinpath("..", "load-data-utils.jl"));
include("plot-utils.jl");

workDir = joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip");

targetDir = joinpath("plots-replicated-with-julia", "calculate_weights");
mkpath(targetDir)

# 1. output weights 
weights = NetCDF.ncread(joinpath(workDir, "weights.nc"), "weight");
models = NetCDF.ncread(joinpath(workDir, "weights.nc"), "model_ensemble");

@assert sum(weights) == 1
@assert length(weights) == 7

xs = 1:length(models)
begin 
    figWeights = Figure();
    ax = Axis(figWeights[1,1], 
        xticks = (xs, models), 
        xticklabelrotation = pi/4,
        xticklabelsize = 10,
        title = "Performance weights",
        xlabel = "Model", 
        ylabel = "Weightws weight (1)"
    );
    barplot!(xs, weights)
    save(joinpath(targetDir, getCurrentTime() * "_performance_weights.png"), figWeights);
end

#### check out how heatmap works ####
# Note: yreversed=true and transposed matrix mat'!
begin
    f = Figure();
    ax = Axis(f[1,1], yreversed=true);
    mat = [0 1 2; 3 0 4; 5 6 0];
    # print(mat)
    hm = heatmap!(ax, mat')
    Colorbar(f[1,2], hm)
    f  
end
#### --------------------------- ####

# Get processed data from work directory that was used for diagnostics
data, models, modelRefs = loadIndependentDistMatrix(workDir, "overall_mean", "overall_mean");
figMean = plotDistMatrices(data, "overall_mean", models, modelRefs)
save(joinpath(targetDir, getCurrentTime() * "_independence_overall_mean.png"), figMean);

data, models, modelRefs = loadIndependentDistMatrix(workDir, "pr_CLIM", "dpr_CLIM");
figPr = plotDistMatrices(data, "pr_CLIM", models, modelRefs)
save(joinpath(targetDir, getCurrentTime() * "_independence_pr_CLIM.png"), figPr);

data, models, modelRefs = loadIndependentDistMatrix(workDir, "tas_CLIM", "dtas_CLIM");
figTas = plotDistMatrices(data, "tas_CLIM", models, modelRefs)
save(joinpath(targetDir, getCurrentTime() * "_independence_tas_CLIM.png"), figTas);



# Root Mean Squared Error plots 
data, models = loadPerformanceDistMatrix(workDir, "overall_mean", "overall_mean");
figPerformOverallMean = plotPerformanceMetric(data, "overall_mean", models)
save(joinpath(targetDir, getCurrentTime() * "_performance_overall_mean.png"), figPerformOverallMean);


data, models = loadPerformanceDistMatrix(workDir, "tas_CLIM", "dtas_CLIM");
figPerformTas = plotPerformanceMetric(data, "tas_CLIM", models)
save(joinpath(targetDir, getCurrentTime() * "_performance_tas_CLIM.png"), figPerformTas);


data, models = loadPerformanceDistMatrix(workDir, "psl_CLIM", "dpsl_CLIM");
figPerformPsl = plotPerformanceMetric(data, "psl_CLIM", models)
save(joinpath(targetDir, getCurrentTime() * "_performance_psl_CLIM.png"), figPerformPsl);


data, models = loadPerformanceDistMatrix(workDir, "pr_CLIM", "dpr_CLIM");
figPerformPr = plotPerformanceMetric(data, "pr_CLIM", models)
save(joinpath(targetDir, getCurrentTime() * "_performance_pr_CLIM.png"), figPerformPr);
