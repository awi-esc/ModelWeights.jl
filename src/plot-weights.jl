"""
    plotPerformanceWeights(wP::DimArray, wP_combined::DimArray)

# Arguments:
- wP: performanceWeights with dimensions: model, variable
- wP_combined: performanceWeights across variables combined into single weight 
per model
"""
function plotPerformanceWeights(
    wP::DimArray, wP_combined::Union{DimArray, Nothing}=nothing, isBarPlot::Bool=true
)
    figures = [];
    models = Array(dims(wP, :model));
    variables = dims(wP, :variable);
    xs = 1:length(models);
    if isBarPlot
        for var in variables
            fig=Figure()
            ax = Axis(fig[1,1], 
                xticks = (xs, models), 
                xticklabelrotation = pi/4,
                xlabel = "Models", 
                ylabel = "Performance weight",
                title = "$var"
            );
            barplot!(ax, xs, Array(wP[variable = At(var)]))
            push!(figures, fig)
        end
    else 
        fig=Figure()
        ax = Axis(fig[1,1], 
            xticks = (xs, models), 
            xticklabelrotation = pi/4,
            xlabel = "Models", 
            ylabel = "Performance weight",
        );
        for var in variables
            scatter!(ax, xs, Array(wP[variable = At(var)]))
            lines!(ax, xs, Array(wP[variable = At(var)]), label = "$var")
        end
        if !isnothing(wP_combined)
            scatter!(ax, xs, Array(wP_combined))
            lines!(ax, xs, Array(wP_combined), label = "combined weights all vars")
        end
        axislegend(ax) 
        push!(figures, fig)
    end
    return figures
end


function plotIndependenceWeights(distances::DimArray)
    figures = []
    models = distances.metadata["full_model_names"]
    if hasdim(distances, :variable)
        for var in dims(distances, :variable)
            fig = plotDistMatrices(
                Array(distances[variable = At(var)]), 
                var, 
                models,
                models
            )
            push!(figures, fig)
        end
    else
        fig = plotDistMatrices(Array(distances), "across variables", models, models)
        push!(figures, fig)
    end
    return figures
end

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
# several figures, one per variable
    # scatter!(
    #     ax, 
    #     Array(independence[variable = At(var)]), 
    #     Array(performance[variable = At(var)]),
    #     label = "$var"
    # )
    # text!(
    #     independence[variable = At(var)], 
    #     performance[variable = At(var)], 
    #     text = Array(dims(performance, :model)), 
    #     align = (:center, :top)
    # )


"""
    plotWeightContributions(independence::DimArray, performance::DimArray)

Make plot of independence vs. performance part of overall weight that combines 
both weights.

The overall weight is computed by dividing the independence part by the 
performance part.

# Arguments: 
- independence: (dims:model) numerator exp(-(D_i/sigma_D)^2)
- performance: (dims:model) denominator 1 + sum_j≠i^M exp(-(S_ij/sigma_S)^2) 
"""
function plotWeightContributions(independence::DimArray, performance::DimArray)
    fig = Figure(size=(800,600), fontsize=16)
    ax = Axis(
        fig[1,1], 
        xlabel =  L"Performance\n $e^{-(D_i/\sigma_D)^2}$", 
        ylabel = L"Independence\n\n $1 + \sum_{j≠i}^M e^{-(S_{ij}/\sigma_S)^2}$", 
        title="(Unnormalized) Performance vs. independence contributions to overall weights",
        xlabelsize = 24,
        ylabelsize = 24
    )

    scatter!(ax, Array(independence), Array(performance))
    text!(
        Array(independence), 
        Array(performance), 
        text = Array(dims(performance, :model)), 
        align = (:center, :top)
    )
    return fig
end


function plotWeights(weights::DimArray)
    fig =  getFigure((16,8), 18);
    models = Array(dims(weights, :model))
    ax = Axis(fig[1,1], 
              xlabel = "Models", 
              ylabel = "weights", 
              xticks = (collect(1:length(models)), models), 
              xticklabelrotation = pi/2);
    xs = 1:length(models);
    scatter!(ax, xs, Array(weights))
    lines!(ax, xs, Array(weights))
    # add line with value if all weights were equal
    n = length(models)
    lines!(ax, xs, [1/n for _ in range(1, n)])
    return fig
end