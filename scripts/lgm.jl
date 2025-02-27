import ModelWeights as mw
using NCDatasets
using DimensionalData
using Setfield
using Distributions
using Statistics
using CairoMakie

# GET THE DATA
# 1. Model data for LGM experiment
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/lgm-cmip5-cmip6";
lgm_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes; 
    subset =  Dict(
        "statistics" => ["CLIM"], 
        "variables" => ["tas"],
        "projects" => ["CMIP5", "CMIP6"], 
        "aliases" => ["lgm"],
        "level_shared_models" => mw.MEMBER
    ),
    preview = false
);
mw.kelvinToCelsius!(lgm_data)
members_lgm = Array(dims(lgm_data["tas_CLIM_lgm"].data, :member));
models_lgm = unique(lgm_data["tas_CLIM_lgm"].data.metadata["model_names"])

# 2. Model data for historical experiment of models with lgm-experiment
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/historical";
# Across variables, only shared model members should be loaded (level_shared_models set to mw.MEMBER) 
# since we want the exact same simulations for all variables when computing weights
# use same model members as for lgm
historical_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes;
    subset = Dict(
        "statistics" =>["CLIM"], 
        "variables" => ["tas"], 
        "projects" => ["CMIP5", "CMIP6"],
        "aliases" => ["historical"], 
        "timeranges" => ["full"],
        "subdirs" => ["20241118", "20250131"],
        "models" => members_lgm,
        "level_shared_models" => mw.MEMBER
    )
);
mw.kelvinToCelsius!(historical_data)
# Load observational data
base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_obs_historical_tas_tos";
config_path = "/albedo/home/brgrus001/ModelWeights/configs/historical_obs";
# aliases and timeranges don't have to match, all data will be loaded that 
# corresponds either to aliases or to timeranges!
obs_data = mw.loadDataFromESMValToolConfigs(
    base_path, config_path;
    dir_per_var = false,
    is_model_data = false,
    subset = Dict(
        "statistics" => ["CLIM"], "variables" => ["tas"], "timeranges" => ["full"]
    )
);
mw.kelvinToCelsius!(obs_data)

# COMPUTATION WEIGHTS
weights_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/";

obs_avg = obs_data["tas_CLIM_historical"].data[source=1]
f = Figure();
mw.plotValsOnMap!(f, obs_data["tas_CLIM_historical"].data[source=1], "Climatology observations")
save("plots/lgm/climatology-observations.png", f)

# 1. Compute weights based on historical data
target_fn = "weights-historical-performance";
target_path = joinpath(weights_dir, target_fn * ".jld2")
config_weights = mw.ConfigWeights(
    performance = Dict("tas_CLIM" => 1),
    independence = Dict("tas_CLIM" => 1),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    ref_perform_weights = "historical", # alias
    #ref_perform_weights = "full", # timerange
    ref_indep_weights = "historical", # alias
    target_path = target_path
);
model_data_vec = collect(values(historical_data));
obs_data_vec = collect(values(obs_data));
dists_perform = mw.computeModelDataRMSE(
    model_data_vec, obs_data_vec, config_weights, "historical"
);
weights = mw.computeWeights(model_data_vec, dists_perform, config_weights)

f1 = mw.plotWeights(weights; title="Weights based on historical data (tas)")
save("plots/lgm/" * target_fn * ".png", f1)

# 2. Apply weights on historical period
df = deepcopy(historical_data);
mw.averageEnsembleMembers!(df)
weighted_avg_hist = mw.applyWeights(df["tas_CLIM_historical"].data, weights.w);
f2 = Figure;
mw.plotValsOnMap!(f2, weighted_avg_hist, "Weighted avg historical tas\n (based on historical data)")
save("plots/lgm/weighted-avg-historical-tas-based-on-historical.png", f2)

###############################################################################
# WEIGHTS BASED ON LGM-COOLING WITH RESPECT TO TIERNEY DATA
# 0. Get data
path_config = "/albedo/home/brgrus001/ModelWeights/configs/ecs-lgm-cooling.yml";
data = mw.loadDataFromYAML(path_config; subset=Dict("level_shared_models" => mw.MEMBER));
# lgm-cooling: here models need to be in same unit for both experiments
mw.kelvinToCelsius!(data)
lgm_cooling = data["tas_CLIM_lgm"].data .- data["tas_CLIM_piControl"].data;
# global lgm-cooling values for each model
global_means = mw.getGlobalMeans(lgm_cooling)

# values from Tierney et al. (2020)
mu_tierney = -6.1
sigma_tierney = 0.4
distr_tierney = Distributions.Normal(mu_tierney, sigma_tierney)
CI_lower(n) = mu_tierney - n * sigma_tierney;
CI_upper(n) = mu_tierney + n * sigma_tierney; 
# 68% CI
l1 = CI_lower(1)
u1 = CI_upper(1)
# 95% CI
l2 = CI_lower(2)
u2 = CI_upper(2)

# plot data
begin
    predictions_gm_tas = Array(global_means);
    ys = repeat([0], length(global_means));
    f3 = Figure(size=(1000,400))
    t="GMST (LGM-piControl)"
    ax = Axis(f3[1,1], title = L"$\Delta$ %$(t)")
    scatter!(ax, predictions_gm_tas, ys, color = :red)

    samples = rand(distr_tierney, 1000);
    hist!(samples; color=(:blue, 0.5), label="Reconstructed Distribution (Data-assimilated) Tierney")
    vlines!(mu_tierney; color=:grey, linewidth=5, label="Mean")
    lines!([l1, u1], [0, 0]; linewidth=5, color=:red, label="68% CI")
    lines!([l2, u2], [0, 0]; linewidth=5, color=:green, alpha=0.5, label="95% CI")

    labels = hasdim(global_means, :model) ? Array(dims(global_means,:model)) : Array(dims(global_means,:member))
    text!(ax, predictions_gm_tas, ys.+5, text=labels, rotation=pi/2)
    axislegend(ax, merge=true, position=:rt)
    f3
end
save("plots/lgm/reconstructed_data_assim_tierney.png", f3)

# 1. compute the weights wrt to lgm data
target_fn = "weights-lgm-cooling";
config_weights = mw.ConfigWeights(
    performance = Dict("tas_lgm-cooling" => 1),
    independence = Dict("tas_CLIM" => 1),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    ref_perform_weights = "lgm",
    ref_indep_weights = "lgm",
    target_path = joinpath(weights_dir, target_fn * ".jld2")
);

dists_perform = mw.getModelLikelihoods(global_means, distr_tierney, "lgm-cooling");
weights = mw.computeWeights([data["tas_CLIM_lgm"]], dists_perform, config_weights)
f4 = mw.plotWeights(weights; title="Weights based on lgm-cooling data wrt Tierney DA-reconstruction (tas)")
save("plots/lgm/weights-lgm-cooling.png", f4)

# these two are identical: (by definition)
mw.plotDistances(weights.performance_distances, "Performance distances";is_bar_plot=false)[1]
mw.plotDistances(dists_perform, "Likelihoods";is_bar_plot=false)[1]

mw.plotDistances(weights.wP, "performance weights"; is_bar_plot = false)[1]


# 2. apply weights on estimated climatological average of historical period
df = deepcopy(data);
mw.averageEnsembleMembers!(df)
weighted_avg_lgm = mw.applyWeights(df["tas_CLIM_historical"].data, weights.w);
f5 = Figure();
mw.plotValsOnMap!(f5, weighted_avg_lgm, "Weighted avg historical tas\n(based on lgm-cooling)")
save("plots/lgm/weighted-avg-historical-tas-based-on-lgm-cooling.png", f5)

# plot difference
f6 = Figure();
mw.plotValsOnMap!(f6, weighted_avg_hist .- weighted_avg_lgm, "tas historical\n Weighted avg based on hist minus weighted avg based on lgm")
save("plots/lgm/weighted-avg-diff.png", f6)

# unweighted average
unweighted_avg = mw.computeWeightedAvg(df["tas_CLIM_historical"].data; use_members_equal_weights = false);
f7 = Figure();
mw.plotValsOnMap!(f7, unweighted_avg, "tas historical\n Unweighted avg")
save("plots/lgm/unweighted-avg.png", f7)

# diff unweighted and weighted
f8 = Figure();
mw.plotValsOnMap!(f8, unweighted_avg .- weighted_avg_lgm, "tas historical\n Unweighted minus weighted based on lgm")
save("plots/lgm/diff_unweighted-minus-weighted-based-on-lgm.png", f8)
f9 = Figure(); 
mw.plotValsOnMap!(f9, unweighted_avg .- weighted_avg_hist, "tas historical\n Unweighted minus weighted based on historical")
save("plots/lgm/diff_unweighted-minus-weighted-based-on-historical.png", f9)




mask = DimArray(zeros(72,36), (Dim{:lon}, Dim{:lat}))
m1 = mw.areaWeightedRMSE(weighted_avg_hist.data, obs_avg, mask)
m2 = mw.areaWeightedRMSE(weighted_avg_lgm.data, obs_avg, mask)



##############################################################################

uncertainties = mw.getUncertaintyRanges(tas_data, weights.w_members);

f3 = mw.plotTempGraph(
    data_graph, 
    (weighted=weighted_avg, unweighted=unweighted_avg),
    uncertainties,
    "Temperature anomaly relative to 1981-2010";
    ylabel = "Temperature anomaly"
)


# simple average across model members
unweighted_means_members = mw.computeWeightedAvg(model_historical_lgm.data["tas_CLIM_lgm"])

# unweighted avg across models
lgm_tas_data = mw.summarizeEnsembleMembersVector(model_historical_lgm.data["tas_CLIM_lgm"], true)
unweighted_means = mw.computeWeightedAvg(lgm_tas_data)
f10 = Figure();
mw.plotValsOnMap!(f10, unweighted_means, "unweighted average LGM: tas_CLIM")
 
# weighted avg across models
weighted_means = mw.computeWeightedAvg(lgm_tas_data; weights=weights.w)
f11 = Figure();
mw.plotValsOnMap!(f11, weighted_means, "weighted means LGM: tas_CLIM")

# weighted and unweighted means should be different:
Array(weighted_means) .== Array(unweighted_means)











# Some experiments with independence weights: How different are they depending on reference period?

# independence weights wrt historical
# independence weights wrt lgm

# add tests that saved config file is adapted when timestamp is added in case file already exists
# @assert target_jld2 == weights_from_disk.config.target_path

########################### 3. PLOTTING ###########################
weights = mw.loadWeightsFromJLD2(joinpath(weights_dir, "weights-lgm-cooling.jld2"))
# if weights loaded from .nc file instead, you can use the following function to load
# respective weights:
# wP = mw.loadWeightsAsDimArray(ds_weights, "wP");
# TODO: Plot generalized distances

# Plot performance and independence weights together
f = mw.plotWeightContributions(weights.wI, weights.wP)

# ds_Sij = ds_weights["Sij"];
# src_names = dimnames(ds_Sij);
# sources = [Array(ds_Sij[src_names[1]]), Array(ds_Sij[src_names[1]])];
# Sij = DimArray(
#     Array(ds_Sij),
#     (Dim{Symbol(src_names[1])}(sources[1]), Dim{Symbol(src_names[2])}(sources[2])), 
#     metadata = Dict(ds_Sij.attrib)
# );
# # Note for weights.Sij dimension names are 'model1' and 'model2', but fn expects twice 'model'
# figs = mw.plotDistancesIndependence(Sij, "model1");
# figs[1]

# figs_performances = mw.plotDistances(weights.performance_distances, "Performance distances");
# figs_performances[1]
# figs_performances[2]

# TODO: fix plotting for independence distances, dimension names are 'model1' and 'model2', but fn expects twice 'model'
#figs = mw.plotDistancesIndependence(weights.independence_distances, "member1");



