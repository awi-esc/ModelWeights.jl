import SimilarityWeights as sw
using DimensionalData
using CairoMakie
using ColorSchemes

# use test data from climwip_test_basic  #
pathToPreprocDir = joinpath(@__DIR__, 
                            "..", 
                            "reproduce-climwip-figs", 
                            "recipe_climwip_test_basic_data", 
                            "preproc", 
                            "calculate_weights_climwip");
climateVariables = ["tas", "pr", "psl"];
sigmaD = 0.5;
sigmaS = 0.5;
##########################################
var_to_preproc_data = Dict{String, Dict{String, String}}();
var_to_preproc_data["CLIM"] = Dict();
for var in climateVariables
    var_to_preproc_data["CLIM"][var] = joinpath(pathToPreprocDir, var * "_CLIM");
end
modelData = sw.loadPreprocData(var_to_preproc_data, ["CMIP"]);
obsData = sw.loadPreprocData(var_to_preproc_data, ["ERA5"]);

# take equal weights for both, performance and independence metric
weightsVarsPerform = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
weightsVarsIndep = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 

wP = sw.generalizedDistancesPerformance(modelData["CLIM"], obsData["CLIM"], weightsVarsPerform);
Di = sw.reduceGeneralizedDistancesVars(wP);

wI = sw.generalizedDistancesIndependence(modelData["CLIM"], weightsVarsIndep);
Sij = sw.reduceGeneralizedDistancesVars(wI);

performances = sw.performanceParts(Di, sigmaD);
independences = sw.independenceParts(Sij, sigmaS);
weights = performances ./ independences;
normalizedWeights = weights ./ sum(weights)


sw.plotPerformanceWeights(wP)
sw.plotPerformanceWeights(wP, Di)

weights = sw.overallWeights(
    modelData, obsData, sigmaD, sigmaS, weightsVarsPerform, weightsVarsIndep
);
# distributed CCSM4 across all 4 CCSM4-models
weightCCSM4 = weights[4]/4;
result = [Array(weights[1:3]); repeat([weightCCSM4], 4)]
models = [Array(dims(weights, :model)); repeat(["CCSM4"], 3)];
weights = DimArray(result, (Dim{:model}(models)))


# plot weights
begin 
    size_inches = (8, 6)
    size_pt = 72 .* size_inches
    fig=Figure(size= size_pt, fontsize=12)

    models = Array(dims(wP, :model));
    xs = 1:length(models);

    ax = Axis(fig[1,1], 
        xticks = (xs, models), 
        xticklabelrotation = pi/4,
        title = "Performance weights",
        xlabel = "Model", 
        ylabel = "Weight"
    );
    #data = Array(wP)
    for (col, var) in enumerate(dims(wP, :variable))
        scatter!(ax, xs, Array(wP[variable = At(var)]))
        lines!(ax, xs, Array(wP[variable = At(var)]), label = "$var")
    end
    
    # add combined weights, summed over all variables
    scatter!(ax, xs, Array(Di))
    lines!(ax, xs, Array(Di), label = "combined perform. weights all vars")

    # add overall weights (indep + perform)
    scatter!(ax, xs, Array(weights))
    lines!(ax, xs, Array(weights), label = "overall weights (with indep.)")

    Legend(fig[1,2], ax)
    fig
end

# plot independence weights
# begin 
#     size_inches = (8, 5)
#     size_pt = 72 .* size_inches
#     fig=Figure(size= size_pt, fontsize=12)

#     models = Array(dims(wP, :model));
#     xs = 1:length(models)
#     ax = Axis(
#         fig[1,1],
#         xticks = (xs, models), 
#         yticks = (xs, models),
#         xticklabelrotation = pi/4,
#         yticklabelrotation = pi/4,
#         title = "Independence weights across variables",
#         yreversed = true
#         )
#     hm = heatmap!(ax, wI_overall, colormap = ColorSchemes.YlGn_4.colors)
#     Colorbar(fig[1, 2], hm)
#     fig
# end


wP = sw.averageEnsembleVector(wP, false);
wI = sw.averageEnsembleMatrix(wI, false);
performance = exp.(-(wP ./ sigmaD).^2);

# note: (+1 in eq. in paper is from when model is compared to itself since exp(0)=1)
independence = dropdims(sum(exp.(-(wI ./ sigmaS).^2), dims=:model2), dims=:model2);
independence = set(independence, :model1 => :model);


tas_ind = independence[variable = At("tas")];
tas_perform = performance[variable = At("tas")];

begin 
    fig=Figure(fontsize=12)
    ax = Axis(fig[1,1], 
              xlabel = "Independence weight", 
              ylabel = "Performance weight", 
              title="Variable: tas"
              )
    scatter!(ax, Array(tas_ind), Array(tas_perform))
    text!(tas_ind, tas_perform, text = Array(dims(tas_perform, :model)), 
          align = (:center, :top))
    fig
end
