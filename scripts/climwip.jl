import ModelWeights as mw
import ModelWeights.Data as mwd
import ModelWeights.Timeseries as mwt
import ModelWeights.Weights as mww
import ModelWeights.Plots as mwp

using CairoMakie
using DimensionalData
using Statistics
using YAXArrays

# get model data
model_ids = mwd.loadModelsFromCSV("../ModelWeightsPaper/work/data/models-brunner-et-al.csv", "Model")
member_ids = mwd.loadModelsFromCSV(
    "../ModelWeightsPaper/work/data/models-brunner-et-al.csv", "Model"; col_variants = "Variants"
)
model_data =  mw.defineDataMap(
    "../ModelWeightsPaper/work/configs/climwip/climwip.yml"; 
    dtype = "cmip", 
    constraint = Dict(
        "level_shared" => "model",
        "models" => member_ids
    )
)
mwd.apply!(model_data, mwt.filterTimeseries, 2014, 2100; 
    ids = ["tas_CLIM-ann_ssp126", "tas_CLIM-ann_ssp585"]
)
mwd.apply!(model_data, mwt.filterTimeseries, 1980, 2014; 
    ids = ["tas_CLIM-ann_historical", "psl_CLIM-ann_historical"], 
    ids_new = ["tas_CLIM-ann_1980-2014", "psl_CLIM-ann_1980-2014"]
)
mwd.apply!(model_data, mwt.filterTimeseries, 1995, 2014; 
    ids = ["tas_CLIM-ann_historical"], 
    ids_new = ["tas_CLIM-ann_1995-2014"]
)

# get observational data
base_dir = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/obs/recipe_ERA5_20250718_180812/preproc/historical"
data_dirs = ["psl_CLIM-ann", "psl_CLIM", "tas_CLIM-ann", "tas_CLIM"]
obs_data = mw.defineDataMap(
    joinpath.(base_dir, data_dirs),
    data_dirs;
    filename_format = :esmvaltool,
    constraint = Dict(
        "variables" => ["tas", "psl"]
    )
)
obs_data = mwd.apply(obs_data, mwd.setDim, :model, nothing, ["ERA5"])
mwd.apply!(obs_data, mwt.filterTimeseries, 1980, 2014; 
    ids = ["tas_CLIM-ann", "psl_CLIM-ann"], 
    ids_new = ["tas_CLIM-ann_1980-2014", "psl_CLIM-ann_1980-2014"]
)
mwd.apply!(obs_data, mwt.filterTimeseries, 1995, 2014; 
    ids = ["tas_CLIM-ann"], 
    ids_new = ["tas_CLIM-ann_1995-2014"]
)

# Compute diagnostics for computing model weights
# make sure that keys of observational datamap for diagnostics are the same as for
# datamap of models
mwd.renameDict!(
    model_data,
    ["tas_CLIM_historical", "psl_CLIM_historical", "tas_CLIM-ann_historical", "psl_CLIM-ann_historical"],
    ["tas_CLIM", "psl_CLIM", "tas_CLIM-ann", "psl_CLIM-ann"]
)

for (i, dm) in enumerate([model_data, obs_data])
    mwd.apply!(
        dm, mwd.anomaliesGM;
        ids = ["tas_CLIM", "psl_CLIM"], 
        ids_new = ["tas_ANOM-GM", "psl_ANOM-GM"]
    )
    mwd.apply!(
        dm, mwt.linearTrend;
        ids = ["tas_CLIM-ann", "psl_CLIM-ann"], 
        ids_new = ["tas_TREND", "psl_TREND"]
    )
    mwd.apply!(
        dm, mwt.detrend;
        ids = ["tas_CLIM-ann", "psl_CLIM-ann"],
        ids_new = ["tas_CLIM-ann-detrended", "psl_CLIM-ann-detrended"]
    )

    for v in ["tas", "psl"]
        dm[v * "_STD"] = mapslices(
            x -> Statistics.std(x), dm[v * "_CLIM-ann-detrended"], dims=(:time,)
        )
    end
end

# Compute Model Weights
config = mww.ConfigWeights(
    performance = Dict(
        "tas_TREND" => 0.5,
        "tas_ANOM-GM" => 0.125,
        "psl_ANOM-GM" => 0.125,
        "tas_STD" => 0.125,
        "psl_STD" => 0.125
    ),
    independence = Dict(
        "tas_CLIM" => 0.5,
        "psl_CLIM" => 0.5
    ),
    sigma_independence = 0.54,
    sigma_performance = 0.43
);


climwip_weights = mww.climwipWeights(model_data, obs_data, config)
fig_weights = mwp.plotWeights(weights.w; sort_by = "combined-climwip")


projections = mwd.apply(
    model_data, mwd.globalMeans;
    ids = ["tas_CLIM-ann", "tas_CLIM-ann_ssp126", "tas_CLIM-ann_ssp585", "tas_CLIM-ann_1995-2014"], 
    ids_new = ["tas_CLIM-ann-GM", "tas_CLIM-ann_ssp126-GM", "tas_CLIM-ann_ssp585-GM", "tas-GM-1995-2014"]
)
mwd.summarizeMembers!(projections)


temp_change_ref = mean(projections["tas-GM-1995-2014"], dims=:time)[time = 1]
for id in ["tas_CLIM-ann-GM", "tas_CLIM-ann_ssp126-GM", "tas_CLIM-ann_ssp585-GM"]
    mwd.apply!(
        projections, mwd.anomalies, temp_change_ref;
        ids = [id],
        ids_new = [id[1:end-3] * "_ANOM-1995-2014"]
    )
end

ids =  ["tas_CLIM-ann_ANOM-1995-2014", "tas_CLIM-ann_ssp126_ANOM-1995-2014", "tas_CLIM-ann_ssp585_ANOM-1995-2014"]
ids_new = ["tas-historical", "tas-ssp126", "tas-ssp585"]
unweighted_avgs = mwd.apply(
    projections, mww.weightedAvg; 
    ids = ids,
    ids_new = ids_new,
    use_members_equal_weights = false
)
unweighted_unc = mwd.apply(projections, mwd.uncertaintyRanges; ids = ids, ids_new = ids_new)

weighted_avgs = mwd.apply(
    projections, mww.weightedAvg;
    ids = ids,
    ids_new = ids_new,
    weights = climwip_weights.w[weight = At("combined-climwip")]
);
# or apply weights for every set of weights, (just performance, just independence, combined)
df = mwd.apply(
    projections, mww.applyWeights, climwip_weights.w;
    ids = ids
)

weighted_unc = mwd.apply(
    projections, mwd.uncertaintyRanges;
    ids = ids,
    ids_new = ids_new,
    w = climwip_weights.w
)

begin
    f = Figure(); 
    ax = Axis(
        f[1,1],
        title = "Near-Surface Air Temperature\nMean and likely (66%) range", 
        xlabel="Year",
        ylabel = "Temperature change (°C) relative to 1995-2014"
    );
    #ylims!(ax, min_val - 0.5, max_val + 5)
    #label_unc_unw(quantiles) = "Non-weighted quantiles: " * join(quantiles, "-")
    #label_unc_w(quantiles) = "Weighted quantiles: " * join(quantiles, "-")

    unweighted_plots = []
    for (id, col) in zip(ids_new, [:grey, :blue, :red])
        plots = mw.Plots.plotTimeseries!(
            ax, 
            unweighted_avgs[id]; 
            uncertainties = coalesce.(unweighted_unc[id], missing => NaN), 
            label = "Non-weighted mean", 
            color_line = col,
            color_unc = :grey,
            linestyle = :dash
        )
        push!(unweighted_plots, plots...)
    end

    weighted_plots = []
    for (id, col) in zip(ids_new, [:grey, :blue, :red])
        plots = mw.Plots.plotTimeseries!(
            ax, 
            weighted_avgs[id]; 
            uncertainties = coalesce.(weighted_unc[id][weight = At("combined-climwip")], missing => NaN), 
            label = "Weighted mean", 
            color_line = col,
            color_unc = col
        )
        push!(weighted_plots, plots...)
    end


    # Dummy plots (invisible, used only for the legend)
    l1 = lines!(ax, 1:2, [NaN, NaN], color = :black, linestyle = :solid)
    l2 = lines!(ax, 1:2, [NaN, NaN], color = :black, linestyle = :dash)
    l3 = lines!(ax, 1:2, [NaN, NaN], color = :red)
    l4 = lines!(ax, 1:2, [NaN, NaN], color = :blue)

    # Create legend with (label, plot_object) pairs
    # Legend(
    #     f[1, 2], 
    #     [l1, l2, l3, l4], ["Weighted", "Unweighted", "SSP5-8.5", "SSP1-2.6"],
    #     title = "Mean and likely (66%) range"    
    # )
    axislegend(
        ax,
        [l1, l2, l3, l4], ["Weighted", "Unweighted", "SSP5-8.5", "SSP1-2.6"];
        position = :lt,
        title = "Mean and likely (66%) range"    
    )
    f
end

mwp.savePlot(f, "../ModelWeightsPaper/work/plots/climwip.png"; overwrite = true)