import SimilarityWeights
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
##########################################
var_to_preproc_data = Dict{String, String}();
for var in climateVariables
    var_to_preproc_data[var] = joinpath(pathToPreprocDir, var * "_CLIM");
end
modelData = SimilarityWeights.loadPreprocData(var_to_preproc_data, ["CMIP"]);
obsData = SimilarityWeights.loadPreprocData(var_to_preproc_data, ["ERA5"]);

# take equal weights for both, performance and independence metric
weightsVarsPerform = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
weightsVarsIndep = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 

wP = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVarsPerform);
wI = SimilarityWeights.getIndependenceWeights(modelData, weightsVarsIndep);

wI_overall = SimilarityWeights.summarizeWeightsAcrossVars(wI);
wP_overall = SimilarityWeights.summarizeWeightsAcrossVars(wP);

SimilarityWeights.plotPerformanceWeights(wP)
SimilarityWeights.plotPerformanceWeights(wP, wP_overall)


sigmaD = 0.5;
sigmaS = 0.5;
weights = SimilarityWeights.getOverallWeights(modelData,
                                              obsData, 
                                              0.5, 
                                              0.5, 
                                              weightsVarsPerform,
                                              weightsVarsIndep
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
    data = Array(wP)
    variables = dims(wP, :variable)
    for (col, var) in enumerate(variables)
        scatter!(ax, xs, data[:, col])
        lines!(ax, xs, data[:, col], label = "$var")
    end
    
    # add combined weights
    scatter!(ax, xs, Array(weights))
    lines!(ax, xs, Array(weights), label = "combined weights all vars")
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


wP = SimilarityWeights.averageEnsembleVector(wP);
wI = SimilarityWeights.averageEnsembleMatrix(wI);
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
