"""
    plotWeights(weights::YAXArray; title::String="", sort_by::String="")

Plot all weights sorted by weight dimension `sort_by` of `weights` if given, otherwise, 
sorted by first value in weight dimension.

Equal weights line is added to plot too.

# Arguments:
- `weights`:
- `title`: 
- `sort_by`: value of weight dimension according to which data is sorted
"""
function plotWeights(weights::YAXArray; title::String="", sort_by::String="")
    sort_by = isempty(sort_by) ? weights.weight[1] : sort_by
    indices = Array(sortperm(weights[weight = At(sort_by)], rev=true))
    
    n_models = length(weights.model)
    xs = 1:n_models
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        #xticks = (xs, Array(weights.model[indices])),
        #xticklabelrotation = pi / 4,
        xlabel = "MODEL",
        ylabel = "Weight",
        title = title,
    )
    labels = Array(weights.weight)
    for label in labels
        ys = Array(weights[weight = At(label)])
        sorted_ys = ys[indices]
        if label == sort_by
            scatterlines!(ax, xs, sorted_ys, label = label, alpha = 0.5)
        else
            scatter!(ax, xs, sorted_ys, label = label, alpha = 0.5)
        end
    end
    # add equal weight for reference
    ys_equal = [1 / n_models for _ in range(1, n_models)]
    lines!(ax, xs, ys_equal, color = :gray, label = "equal weighting", linestyle = :dashdot)
    axislegend(ax, position = :rt, merge = true)
    return fig
end



function plotDistancesByVar(dists::AbstractArray, title::String; is_bar_plot::Bool = true)
    dists_var = reduce(+, dists, dims = :diagnostic)
    dists_var = DimensionalData.set(
        dists_var,
        :diagnostic => ["sum of diagnostics: $(Array(dims(dists, :diagnostic))...)"],
    )
    plotDistances(dists_var, title; is_bar_plot = is_bar_plot)
end


"""
    plotDistances(dists::AbstractArray, title::String; is_bar_plot::Bool=true)

Plot figure of distances for every combination of variable and diagnostic in 
`dists`.

# Arguments:
- `dists`: must have dimensions `:variable`, `:diagnostic`
- `title`:
- `is_bar_plot`:
"""
function plotDistances(dists::AbstractArray, title::String; is_bar_plot::Bool = true)
    models = hasdim(dists, :member) ? dims(dists, :member) : dims(dists, :model)
    variables = dims(dists, :variable)
    # if isnothing(variables)
    #     variables = ["variables combined"]
    #     dists_reshaped = reshape(dists, (size(dists)..., 1))
    #     dists = DimArray(dists_reshaped, (dims(dists)..., Dim{:variable}(variables)))
    # end
    diagnostics = dims(dists, :diagnostic)
    # if isnothing(diagnostics)
    #     diagnostics = ["Generalized distance"]
    #     dists_reshaped = reshape(dists, (size(dists)..., 1))
    #     dists = DimArray(dists_reshaped, (dims(dists)..., Dim{:diagnostic}(Array(diagnostics))))
    # end

    figures = []
    xs = 1:length(models)
    for diag in Array(diagnostics)
        for var in variables
            fig = Figure()
            ax = Axis(
                fig[1, 1],
                xticks = (xs, Array(models)),
                xticklabelrotation = pi / 4,
                xlabel = "Model member",
                title = title * "Variable: $var, Diagnostic: $diag",
            )
            ys = vec(dists[variable=At(var), diagnostic=At(diag)])
            if is_bar_plot
                barplot!(ax, xs, ys, label = "$var")
            else
                scatter!(ax, xs, ys)
                lines!(ax, xs, ys, label = "$var")
            end
            axislegend(ax, merge = true, position = :lt)
            push!(figures, fig)
        end
    end
    return figures
end


"""
    plotDistancesIndependence(distances::AbstractArray, dimname::String)

# Arguments:
- `distances`:
- `dimname`:
"""
function plotDistancesIndependence(distances::AbstractArray, dimname::String)
    figures = []
    ensembles = collect(dims(distances, Symbol(dimname)))
    if hasdim(distances, :variable)
        for var in dims(distances, :variable)
            for diag in dims(distances, :diagnostic)
                fig = plotDistMatrices(
                    Array(distances[variable=At(var), diagnostic=At(diag)]),
                    var * "_" * diag,
                    ensembles,
                    ensembles,
                )
                push!(figures, fig)
            end
        end
    else
        fig = plotDistMatrices(
            Array(distances),
            "Generalized distances Sij",
            ensembles,
            ensembles,
        )
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
    xlabel::String = "Normalized performance weights",
    ylabel::String = "Normalized independence weights",
    title::String = "",
)
    #L"Performance\n $e^{-(D_i/\sigma_D)^2}$",
    #L"Independence\n\n $1 + \sum_{j≠i}^M e^{-(S_{ij}/\sigma_S)^2}$",  
    fig = Figure(size = (800, 600), fontsize = 16)
    ax = Axis(
        fig[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xlabelsize = 24,
        ylabelsize = 24,
    )
    scatter!(ax, Array(performance), Array(independence))
    m = maximum([maximum(independence), maximum(performance)])
    extra = 0.0005
    lines!(ax, [0, m + extra], [0, m + extra], color = :gray)
    text!(
        Array(performance),
        Array(independence),
        text = Array(dims(performance, :model)),
        align = (:center, :top),
    )
    return fig
end
