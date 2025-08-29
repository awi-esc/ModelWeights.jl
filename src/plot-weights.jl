"""
    plotWeights(weights::YAXArray; title::String="", sort_by::String="")

Plot `weights` sorted by weight dimension `sort_by` of `weights` if given, otherwise, sorted by first value in weight dimension.

Equal weights line is added to the plot.

# Arguments:
- `weights::YAXArray`: with dimensions :weight and :model
- `title::String`:
- `sort_by::String`: value of weight dimension according to which data is sorted
"""
function plotWeights(
    weights::YAXArray; 
    title::String = "", 
    sort_by::String = "", 
    fs::Number = 15, 
    legend_inside::Bool = true,
    leg_pos::Symbol = :rt,
    leg_frame::Bool = true,
    leg_orientation::Symbol = :vertical,
    leg_rows::Int = 1
)
    sort_by = isempty(sort_by) ? weights.weight[1] : sort_by
    indices = Array(sortperm(weights[weight = At(sort_by)], rev=true))
    
    n_models = length(weights.model)
    xs = 1:n_models
    fig = Figure(size=(600,350))
    ax = Axis(
        fig[1, 1],
        xticks = (xs, Array(weights.model[indices])),
        xticklabelrotation = pi/2,
        xlabel = "Model",
        ylabel = "Weight",
        title = title,
        titlesize = fs,
        xticklabelsize = fs,
        yticklabelsize = fs,
        xlabelsize = fs,
        ylabelsize = fs
    )
    labels = lookup(weights, :weight)
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
    lines!(ax, xs, ys_equal, color = :gray, label = "equal weighting", linestyle=:dash)
    if leg_orientation == :horizontal
        leg = axislegend(
            ax, 
            position = leg_pos,
            merge = true, 
            framevisible = leg_frame, 
            orientation = leg_orientation,
            nbanks = leg_rows
        )
    else
        leg = axislegend(
            ax, 
            position = leg_pos,
            merge = true, 
            framevisible = leg_frame, 
            orientation = leg_orientation
        )
    end
    if !legend_inside
        fig[1,2] = leg
    end 
    return fig
end


"""
    plotDistances(dists::AbstractArray, title::String; is_bar_plot::Bool=true)

Plot figure of distances.

# Arguments:
- `dists::AbstractArray`:
- `title::String`:
- `is_bar_plot::Bool = true`:
"""
function plotDistances(dists::AbstractArray, title::String; is_bar_plot::Bool = true)
    model_dim = Data.modelDim(dists)
    models = lookup(dists, model_dim)
    xs = 1:length(models)
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xticks = (xs, collect(models)),
        xticklabelrotation = pi / 4,
        xlabel = String(model_dim),
        title = title
    )
    ys = collect(dists.data)
    if is_bar_plot
        barplot!(ax, xs, ys)
    else
        scatter!(ax, xs, ys)
        lines!(ax, xs, ys)
    end
    return fig
end


"""
    plotDistancesIndependence(
        distances::AbstractArray, dimname::String; title::String="Generalized distances Sij"
    )

# Arguments:
- `distances::AbstractArray`:
- `dimname::String`:
- `title::String`:
"""
function plotDistancesIndependence(
    distances::AbstractArray, dimname::String; title::String="Generalized distances Sij"
)
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
            title,
            ensembles,
            ensembles,
        )
        push!(figures, fig)
    end
    return figures
end
