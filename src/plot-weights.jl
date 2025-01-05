"""
    plotWeights(
    w::DimArray; 
    ylabel::String="weight", is_bar_plot::Bool=true, dimname::String="model"
)

# Arguments:
- `w`: performanceWeights with dimensions: model/member, variable
- `ylabel`: y-axis label
- `is_bar_plot`: if false scatter plot returned
- `dimname`: dimension name
"""
function plotWeights(
    w::DimArray; 
    ylabel::String="weight", is_bar_plot::Bool=true, dimname::String="model"
)
    figures = [];
    models = Array(dims(w, Symbol(dimname)))
    n_models = length(models)
    variables = dims(w, :variable)
    hasDimVariable = true
    if isnothing(variables)
        hasDimVariable = false
        variables = ["variables combined"]
        w_reshaped = reshape(w, (size(w)..., 1))
        w = DimArray(w_reshaped, (dims(w)..., Dim{:variable}(variables)))
    end
    xs = 1 : n_models
    if is_bar_plot
        for var in variables
            fig=Figure()
            ax = Axis(fig[1,1], 
                xticks = (xs, models), 
                xticklabelrotation = pi/4,
                xlabel = dimname,
                ylabel = ylabel,
                title = "$var"
            );
            barplot!(ax, xs, vec(w[variable = At(var)]))
            push!(figures, fig)
        end
    else 
        fig=Figure()
        ax = Axis(fig[1,1], 
            xticks = (xs, models), 
            xticklabelrotation = pi/4,
            xlabel = uppercase(dimname), 
            ylabel = ylabel,
        );
        for var in variables
            scatter!(ax, xs, Array(w[variable = At(var)]))
            lines!(ax, xs, Array(w[variable = At(var)]), label = "$var")
        end
        if hasDimVariable
            axislegend(ax)
        else
            # add equal weight for reference
            lines!(ax, xs, [1/n_models for _ in range(1, n_models)])
        end 
        push!(figures, fig)
    end
    return figures
end


"""
    plotDistancesPerformance(
    dists::DimArray; is_bar_plot::Bool=true, dimname::String="model"
)
"""
function plotDistancesPerformance(dists::DimArray; is_bar_plot::Bool=true)
    dimname = hasdim(dists, :model) ? :model : :member
    figs = plotWeights(
        dists; ylabel = "Distances performance",
        is_bar_plot = is_bar_plot, 
        dimname = String(dimname)
    )
    return figs
end

function plotDistancesIndependence(distances::DimArray; dimname::String="model")
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
- independence: (dims:model) normalized independence weights
- performance: (dims:model) normalized performance weights 
"""
function plotWeightContributions(
    independence::DimArray,
    performance::DimArray;
    xlabel::String="Normalized performance weights", 
    ylabel::String="Normalized independence weights",
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


# """
#     plotWeights(weights::DimArray; title::String="")

# # Arguments:
# - `weights`: overall weights, one for each ensemble
# - `title`: plot title
# """
# function plotWeights(weights::DimArray; title::String="")
#     fig =  getFigure((16,8), 18);
#     models = Array(dims(weights, :model))
#     ax = Axis(fig[1,1], 
#               xlabel = "Models", 
#               ylabel = "weights",
#               title = title,
#               xticks = (collect(1:length(models)), models), 
#               xticklabelrotation = pi/2,
#     );
#     ylims!(ax, 0, 1)
#     xs = 1:length(models);
#     scatter!(ax, xs, Array(weights))
#     lines!(ax, xs, Array(weights))
#     # add line with value if all weights were equal
#     n = length(models)
#     lines!(ax, xs, [1/n for _ in range(1, n)])
#     return fig
# end