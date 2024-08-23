"""
    plotPerformanceWeights(wP::DimArray)

# Arguments:
- wP: performanceWeights with dimensions: model, variable
- wP_combined: performanceWeights across variables combined into single weight 
per model
"""
function plotPerformanceWeights(wP::DimArray, wP_combined::DimArray)
    fig=Figure()
    models = Array(dims(wP, :model));
    xs = 1:length(models);
    ax = Axis(fig[1,1], 
        xticks = (xs, models), 
        xticklabelrotation = pi/4,
        xlabel = "Models", 
        ylabel = "Performance weight"
    );
    data = Array(wP)
    variables = dims(wP, :variable)
    for (col, var) in enumerate(variables)
        scatter!(ax, xs, data[:, col])
        lines!(ax, xs, data[:, col], label = "$var")
    end
    scatter!(ax, xs, Array(wP_combined))
    lines!(ax, xs, Array(wP_combined), label = "combined weights all vars")
    axislegend(ax)
    return fig
end

"""
    plotWeightContributions(independence::DimArray, performance::DimArray)

Make plot of independence vs. performance part of overall weight that combines 
both weights.

The overall weight is computed by dividing the independence part by the 
performance part.

# Arguments: 
- independence: (dims:model x variable) numerator exp(-(D_i/sigma_D)^2)
- performance: (dims:model x variable) denominator 1 + sum_j≠i^M exp(-(S_ij/sigma_S)^2) 
"""
function plotWeightContributions(independence::DimArray, performance::DimArray)
    fig = Figure(size=(800,600), fontsize=16)
    ax = Axis(fig[1,1], 
              xlabel =  L"Independence\n $e^{-(D_i/\sigma_D)^2}$", 
              ylabel = L"Performance\n\n $1 + \sum_{j≠i}^M e^{-(S_{ij}/\sigma_S)^2}$", 
              title="Independence vs. performance contributions to overall weights",
              xlabelsize = 24,
              ylabelsize = 24
              )
    variables = dims(independence, :variable)
    for (col, var) in enumerate(variables)
        scatter!(ax, Array(independence)[:, col], Array(performance)[:, col],
                 label = "$var")
        text!(independence[:, col], performance[:, col], 
              text = Array(dims(performance, :model)), 
              align = (:center, :top))
    end
    axislegend(ax, position = :ct, orientation = :horizontal)
    return fig
end