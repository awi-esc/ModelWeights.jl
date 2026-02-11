"""
    plotWeights(
        weights::YAXArray; 
        fs::Number = 15,
        sort_by::String = "",
        weights_labels::Union{Dict{String, String}, Nothing} = nothing, 
        fig_size::Tuple = (600,450)
        sort_by::String=""
    )

Plot `weights` sorted by weight dimension `sort_by` of `weights` if given, otherwise, sorted by first value in weight dimension.

Equal weights line is added to the plot.

# Arguments:
- `weights::YAXArray`: with dimensions :weight and :model
- `fs::Number`: fontsize
- `sort_by::String`: value of weight dimension according to which data is sorted
- `weights_labels::Union{Dict{String, String}, Nothing} = nothing`: titles of subplots. Default is dimension names for dimension :weight
- `fig_size::Tuple`:
"""
function plotWeights(
    weights::YAXArray; 
    fs::Number = 15,
    sort_by::String = "",
    labels::Union{AbstractArray{String}, Nothing} = nothing, 
    fig_size::Tuple = (600,450)
)
    if !hasdim(weights, :weight)
        indices = 1:length(weights)
    else
        sort_by = isempty(sort_by) ? weights.weight[1] : sort_by
        indices = Array(sortperm(weights[weight = At(sort_by)], rev=true))
    end

    n_models = length(weights.model)
    xs = 1:n_models
    fig = Figure(size=fig_size)
    ys_equal = [1 / n_models for _ in range(1, n_models)]

    weights_labels = lookup(weights, :weight) 
    if isnothing(labels)
        labels = weights_labels
    end

    # axis for last row has model names in xaxis labels
    n_weights = length(weights_labels)
    ax = Axis(
        fig[n_weights, 1],
        xticks = (xs, Array(weights.model)[indices]),
        xticklabelrotation = pi/2,
        xlabel = "Model",
        title = labels[n_weights],
        titlesize = fs,
        xticklabelsize = fs-2,
        yticklabelsize = fs-2,
        xlabelsize = fs,
        ylabelsize = fs
    )
    
    for (i, label) in enumerate(weights_labels)
        ax_i = i==n_weights ? ax : Axis(fig[i,1], title=labels[i])
        ys = Array(weights[weight = At(label)])
        sorted_ys = ys[indices]
        scatterlines!(ax_i, xs, sorted_ys, alpha = 0.5)
        lines!(ax_i, xs, ys_equal, color = :gray, label = "equal weighting", linestyle=:dash)
    end
    
    Label(fig[:,0], "Weight", rotation=pi/2)

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
    ylabel::String = "CRPSS weighted vs. unweighted",
    leg_pos::Symbol = :rb,
    fs::Number = 12
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
                xticklabelrotation = pi / 2,
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
                axislegend(
                    ax, 
                    [PolyElement(color=c) for c in colors], 
                    labels, 
                    position=leg_pos, 
                    merge=true
                )
            end
            push!(figs, fig)
        end
    end
    return figs
end


function crpssBoxPlot(
    crpss_wIP::YAXArray,  colors::Vector{Symbol}; 
    title::String="", leg_pos::Symbol=:rb, ylabel::String="CRPSS", 
    yticks::AbstractArray{<:Number},
    ymax::Number=1
)
    tps = lookup(crpss_wIP, :target_period)
    scenarios = lookup(crpss_wIP, :scenario)
    f = Figure()
    ax = Axis(
        f[1, 1],
        xticks = (1:length(tps), tps),
        yticks = (yticks, string.(yticks)),
        ylabel = ylabel,
        title = title
    )
    ylims!(ax, nothing, ymax)
    hlines!(ax, 0; color=:black)
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
        merge=true,
        framevisible=false
    )
    return f
end


"""
    boxplotMCMCWeights(
        chains, 
        models;
        chain = 1,
        xticks::AbstractArray=0.1:0.1:1,
        xlims::Union{Tuple, Nothing}=nothing,
        title::String="",
        fig_size=(600,400)
    )

# Arguments:
- `chains`: a vector where each element corresponds to one MCMC chain, which contain 
matrices of size N_iter x N_models
- `models`: model names
"""
function boxplotMCMCWeights(
    chains, 
    models;
    chain = 1,
    xticks::AbstractArray=0.1:0.1:1,
    xlims::Union{Tuple, Nothing}=nothing,
    title::String="",
    fig_size=(600,400)
)
    N_models = length(models)
    samples = chains[chain]
    N_iter = size(samples, 1)
    f = Figure(size=fig_size)
    offset = 1
    ys = [(x-1) * offset for x in 1:N_models]
    ax = Axis(f[1,1], xticks = (xticks, string.(xticks)), yticks = (ys, models), title=title, xlabel="Weight")
    Makie.ylims!(ax, -offset/2, maximum(ys) + offset)
    if !isnothing(xlims)
        Makie.xlims!(ax, xlims...)
    end
    for m in 1:length(models)
        Makie.boxplot!(
            ax, 
            fill(ys[m], N_iter), samples[:, m], 
            orientation = :horizontal,
            color=:grey
        )
        Makie.scatter!(
            ax, 
            mean(samples[:, m]),
            ys[m],
            color=RGBf(206/255, 250/255, 220/255),
            marker = :xcross,
            markersize=14,
            label = "BMA-mean"
        )
    end
    Makie.vlines!(ax, 1/N_models, color=:grey, linestyle=:dash)
    return f
end


"""
    plotCorrWeights(models_labels, weights)

# Arguments:
- `models_labels::AbstractVector`: names of models
- `weights::AbstractArray`: samples (iterations) x models
- `pairs::AbstractVector`: tuples for models to be compared against one another as indices
"""
function plotCorrWeights(
    models_labels::AbstractVector,
    weights::AbstractArray, 
    pairs::AbstractVector;
    color::Symbol = :grey,
    color_means::Union{Symbol, ColorTypes.RGB, Nothing} = CairoMakie.RGBf(206/255, 250/255, 220/255),
    fig_size::Union{Tuple, Nothing} = nothing
)
    f = isnothing(fig_size) ? Figure() : Figure(size=fig_size)
    for (i, (m1, m2)) in enumerate(pairs)
        ax = Axis(
            f[1, i], xlabel = models_labels[m1], ylabel = models_labels[m2],
            xticks = (0:0.1:1, string.(0:0.1:1)),
            yticks = (0:0.1:1, string.(0:0.1:1))
        );
        Makie.scatter!(ax, weights[:, m1], weights[:, m2], color=color)
        if !isnothing(color_means)
            Makie.scatter!(ax, mean(weights[:, m1]), mean(weights[:, m2]), color=color_means, marker='x', markersize=30)
        end
    end
    return f
end
