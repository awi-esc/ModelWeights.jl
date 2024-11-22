"""
    plotPerformanceWeights(
    wP::DimArray; label::String="Performance Weight", isBarPlot::Bool=true
)

# Arguments:
- `wP`: performanceWeights with dimensions: ensemble, variable
- `isBarPlot`: if true barplot, else scatter plot returned
"""
function plotPerformanceWeights(
    wP::DimArray; label::String="Performance Weight", isBarPlot::Bool=true, dimname::String="model"
)
    figures = [];
    models = Array(dims(wP, Symbol(dimname)));
    variables = dims(wP, :variable);
    hasDimVariable = true
    if isnothing(variables)
        hasDimVariable = false
        variables = ["variables combined"]
        wP_reshaped = reshape(wP, (size(wP)..., 1))
        wP = DimArray(wP_reshaped, (dims(wP)..., Dim{:variable}(variables)))
    end
    xs = 1:length(models);
    if isBarPlot
        for var in variables
            fig=Figure()
            ax = Axis(fig[1,1], 
                xticks = (xs, models), 
                xticklabelrotation = pi/4,
                xlabel = uppercase(dimname),
                ylabel = label,
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
            xlabel = uppercase(dimname), 
            ylabel = label,
        );
        for var in variables
            scatter!(ax, xs, Array(wP[variable = At(var)]))
            lines!(ax, xs, Array(wP[variable = At(var)]), label = "$var")
        end
        if hasDimVariable
            axislegend(ax)
        end 
        push!(figures, fig)
    end
    return figures
end


function plotIndependenceWeights(distances::DimArray; dimname::String="model")
    figures = []
    ensembles = collect(dims(distances, Symbol(dimname)))
    if hasdim(distances, :variable)
        for var in dims(distances, :variable)
            fig = plotDistMatrices(
                Array(distances[variable = At(var)]), 
                var, 
                ensembles,
                ensembles
            )
            push!(figures, fig)
        end
    else
        fig = plotDistMatrices(Array(distances), "across variables", ensembles, ensembles)
        push!(figures, fig)
    end
    return figures
end


"""
    plotWeightContributions(independence::DimArray, performance::DimArray)

Plot performance against independence weights.

# Arguments: 
- independence: (dims:ensemble) numerator exp(-(D_i/sigma_D)^2)
- performance: (dims:ensemble) denominator 1 + sum_j≠i^M exp(-(S_ij/sigma_S)^2) 
"""
function plotWeightContributions(
    independence::DimArray,
    performance::DimArray;
    xlabel::String="Performance", 
    ylabel::String="Independence",
    title::String=""
)
    #L"Performance\n $e^{-(D_i/\sigma_D)^2}$",
    #L"Independence\n\n $1 + \sum_{j≠i}^M e^{-(S_{ij}/\sigma_S)^2}$",  
    fig = Figure(size=(800,600), fontsize=16)
    ax = Axis(
        fig[1,1], xlabel = xlabel, ylabel = ylabel, title=title, 
        xlabelsize = 24, ylabelsize = 24
    )
    scatter!(ax, Array(performance), Array(independence))
    m = maximum([maximum(independence), maximum(performance)])
    extra = 0.0005
    lines!(ax, [0, m+extra], [0, m+extra], color=:gray)
    text!(
        Array(performance), 
        Array(independence), 
        text = Array(dims(performance, :model)), 
        align = (:center, :top)
    )
    return fig
end


"""
    plotWeights(weights::DimArray; title::String="")

# Arguments:
- `weights`: overall weights, one for each ensemble
- `title`: plot title
"""
function plotWeights(weights::DimArray; title::String="")
    fig =  getFigure((16,8), 18);
    models = Array(dims(weights, :model))
    ax = Axis(fig[1,1], 
              xlabel = "Models", 
              ylabel = "weights",
              title = title,
              xticks = (collect(1:length(models)), models), 
              xticklabelrotation = pi/2,
    );
    ylims!(ax, 0, 1)
    xs = 1:length(models);
    scatter!(ax, xs, Array(weights))
    lines!(ax, xs, Array(weights))
    # add line with value if all weights were equal
    n = length(models)
    lines!(ax, xs, [1/n for _ in range(1, n)])
    return fig
end