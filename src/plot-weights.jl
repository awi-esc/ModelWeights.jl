"""
    plotWeights(weights::Weights; title::String="")

Plot all weights together, i.e. overall weights, performance-only and 
independence-only weights.

# Arguments:
- `weights`:
- `title`: 
"""
function plotWeights(weights::Weights; title::String="")
    models = Array(dims(weights.wP, :model))
    n_models = length(models) 
    xs = 1 : n_models
    fig=Figure()
    ax = Axis(fig[1,1], 
        xticks = (xs, models), 
        xticklabelrotation = pi/4,
        xlabel = "MODEL", 
        ylabel = "Weight",
        title = title
    );
    ys = [weights.wP, weights.wI, weights.w]
    colors = [:red, :green, :blue]
    labels = ["performance", "independence", "overall"]
    for idx in 1:3
        scatter!(ax, xs, Array(ys[idx]), color=colors[idx], label=labels[idx], alpha=0.5)
        lines!(ax, xs, Array(ys[idx]), color=colors[idx], label=labels[idx], alpha=0.5)
    end
    # add equal weight for reference
    lines!(ax, xs, [1/n_models for _ in range(1, n_models)], 
        color=:gray, label="unweighted", linestyle=:dashdot)
    axislegend(ax, position=:lt, merge = true)
    return fig
end


function plotDistances(dists::DimArray, title::String; is_bar_plot::Bool=true)
    models = dims(dists, :member)
    if isnothing(models)
        models = Array(dims(dists, :model))
    end
    variables = dims(dists, :variable)
    if isnothing(variables)
        variables = ["variables combined"]
        dists_reshaped = reshape(dists, (size(dists)..., 1))
        dists = DimArray(dists_reshaped, (dims(dists)..., Dim{:variable}(variables)))
    end
    diagnostics = dims(dists, :diagnostic)
    if isnothing(diagnostics)
        diagnostics = ["Generalized distance"]
        dists_reshaped = reshape(dists, (size(dists)..., 1))
        dists = DimArray(dists_reshaped, (dims(dists)..., Dim{:diagnostic}(Array(diagnostics))))
    end

    figures = []
    xs = 1 : length(models)
    for diag in Array(diagnostics)
        fig=Figure()
        ax = Axis(
            fig[1,1], xticks = (xs, Array(models)), xticklabelrotation = pi/4,
            xlabel = "Model member", title = title * " -- Diagnostic: $diag"
        );
        for var in variables
            ys = vec(dists[variable = At(var), diagnostic = At(diag)])
            if is_bar_plot
                barplot!(ax, xs, ys, label = "$var")
            else
                scatter!(ax, xs, ys)
                lines!(ax, xs, ys, label = "$var")
            end
        end
        axislegend(ax, merge=true, position=:ct)
        push!(figures, fig)
    end
    return figures
end


"""
    plotDistancesIndependence(distances::DimArray, dimname::String)

# Arguments:
- `distances`:
- `dimname`:
"""
function plotDistancesIndependence(distances::DimArray, dimname::String)
    figures = []
    ensembles = collect(dims(distances, Symbol(dimname)))
    if hasdim(distances, :variable)
        for var in dims(distances, :variable)
            for diag in dims(distances, :diagnostic)
                fig = plotDistMatrices(
                    Array(distances[variable = At(var), diagnostic=At(diag)]), 
                    var * "_" * diag, ensembles, ensembles
                )
                push!(figures, fig)
            end
        end
    else
        fig = plotDistMatrices(Array(distances), "Generalized distances Sij", ensembles, ensembles)
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