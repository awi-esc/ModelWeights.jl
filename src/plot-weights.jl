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


function plotCRPSSPseudoObs(
    crpss_all::Vector{T}, labels::Vector{String}, colors::Vector{Symbol};
    ylabel::String = "CRPSS forecast: weighted, baseline: unweighted"
) where {T <: Any}
    figs = []
    models = sort(lookup(crpss_all[1], :pseudo_obs))
    n = length(models)
    target_periods = lookup(crpss_all[1], :target_period)
    scenarios = lookup(crpss_all[1], :scenario)
    for tp in target_periods
        for scenario in scenarios
            # get data for respective scenario and target time period:
            title = scenario * " for target period " * tp
            xs = 1:n
            fig = Figure()
            ax = Axis(
                fig[1, 1],
                xticks = (xs, models),
                xticklabelrotation = pi / 4,
                ylabel = ylabel,
                title = title
            )
            n_groups = length(crpss_all)
            crpss_vals_all = map(
                i -> map(
                        m -> only(crpss_all[i][scenario=At(scenario), target_period=At(tp), pseudo_obs=At(m)]),
                        models
                    ),
                1:n_groups
            )
            groups = vcat(map(i -> repeat([i], length(xs)), 1:n_groups)...)
            xs_flat = vcat(repeat(xs, n_groups)...)
            barplot!(
                ax, xs_flat, vcat(crpss_vals_all...); 
                dodge = groups,
                color = repeat(colors, inner=n),
                label = repeat(labels, inner=n)
            )
            if length(labels) > 1
                axislegend(ax, [PolyElement(color=c) for c in colors], labels, position=:rb, merge=true)
            end
            push!(figs, fig)
            #name = join(["CRPSS-combined-vs-classic-wP", scenario, tp, ".png"], "-", "")
            #mwp.savePlot(fig, joinpath(plot_dir, name); overwrite=true)
        end
    end
    return figs
end


function crpssBoxPlot(crpss_wIP::YAXArray; title::String="", leg_pos::Symbol=:rb)
    tps = lookup(crpss_wIP, :target_period)
    scenarios = lookup(crpss_wIP, :scenario)
    colors = [:sienna2, :orchid3]
    f = Figure()
    ax = Axis(
        f[1, 1],
        xticks = (1:length(tps), tps),
        ylabel = "CRPSS",
        title = title
    )
    # dodged barplot
    ys_all = []
    grps_all = []
    xs_all = []
    for (i, tp) in enumerate(tps)
        ys = collect.(map(sc -> crpss_wIP[target_period = At(tp), scenario=At(sc)], scenarios))
        push!(ys_all, vcat(ys...))
        grps = vcat(map(x -> repeat([x], length(ys[x])), 1:length(scenarios))...)
        push!(grps_all, grps)
        xs = repeat([i], length(vcat(ys...)))
        push!(xs_all, xs)
        
        boxplot!(
            ax, xs, vcat(ys...); 
            dodge = grps,
            color = repeat(colors, inner=length(ys[1])),
            label = repeat(scenarios, inner=length(ys[1]))
        )
    end
    axislegend(
        ax, 
        [PolyElement(color=c) for c in colors], 
        scenarios, 
        position=leg_pos,
        orientation = :horizontal,
        merge=true
    )
    return f
end