import SimilarityWeights as sw
using NCDatasets
using DimensionalData

########################### 1. LOADING DATA ###########################
# Model data, for LGM experiments
base_path =  "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM/";
config_path = "/albedo/home/brgrus001/SimilarityWeights/configs/lgm-cmip5-cmip6";
model_data_lgm = sw.loadData(
    base_path, 
    config_path;
    dir_per_var = true,
    common_models_across_vars = true,
    subset = Dict(
        "statistic" => ["CLIM"],
        "variable" => ["tas", "tos"],
        "projects" => ["CMIP5", "CMIP6"],
        "models" => Vector{String}(), # same as not setting it
        "subdirs" => ["20241114"] # if dir_per_var is true only subdirs containing any are considered
        )
);
#model_data_lgm_shared = sw.getCommonModelsAcrossVars(model_data_lgm, :member);
model_members_lgm = Array(dims(first(values(model_data_lgm.data)), :member))
models_lgm  = unique(first(values(model_data_lgm.data)).metadata["model_names"])


# Model data for historical experiments for those models from above that also have lgm - experiments
base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/";
config_path = "/albedo/home/brgrus001/SimilarityWeights/configs/historical";
model_data_historical = sw.loadData(
    base_path, config_path;
    dir_per_var = true,
    common_models_across_vars = true,
    subset = Dict(
        "statistic" => ["CLIM"],
        "variable" => ["tas", "tos"],
        "alias" => ["historical", "historical0"],
        "timerange" => ["full"],
        "projects" => ["CMIP5", "CMIP6"],
        "models" => model_members_lgm,
        #"models" => models_lgm,
        "subdirs" => ["20241121", "20241118"] # if dir_per_var is true only subdirs containing any are considered
    )
);
#model_data_historical = sw.getCommonModelsAcrossVars(model_data_historical, :member);

# Observational data
obs_data = sw.loadData(
    "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_obs_historical_20241119_124434",
    "/albedo/home/brgrus001/SimilarityWeights/configs/historical_obs";
    dir_per_var = false,
    is_model_data = false,
    subset = Dict(
        "statistic" => ["CLIM"],
        "variable" => ["tas", "tos"],
        "alias" => ["historical", "historical3"],
        "projects" => ["ERA5"] # default value (same as not setting it)
        #"timeranges" => ["1980-2014"]
    )
);


########################### 2. COMPUTATION WEIGHTS ###########################
# if target_dir is provided within config_weights, the weights will directly 
# be saved and written to a file
config_weights = sw.ConfigWeights(
    performance = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    independence = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    ref_period = "historical", # alias
    #ref_period = "full", # timerange
    target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/"
);
weights = sw.computeWeights(model_data_historical, obs_data, config_weights);

# weights can also be  saved separately:
target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/"
target_fn = "weights-lgm-models.nc"
sw.saveWeights(weights, target_dir; target_fn = target_fn)


########################### 3. PLOTTING ###########################
# Plot weights/generalized distances
path_weights = joinpath(target_dir, target_fn);
ds_weights = NCDataset(path_weights);


ds_weights = NCDataset("/albedo/work/projects/p_pool_clim_data/britta/weights/weights_2024-11-27_09_22.nc")
# Plot performance weights
wP = sw.loadWeightsAsDimArray(ds_weights, "wP");
fig_wP, = sw.plotWeights(wP; isBarPlot=false, label="performance weight");
fig_wP

# Plot independence weights
wI = sw.loadWeightsAsDimArray(ds_weights, "wI");
fig_wI, = sw.plotWeights(wI; isBarPlot=false, label="independence weight");
fig_wI

# Plot performance and independence weights together
f = sw.plotWeightContributions(wI, wP)

# Plot overall weights
w = sw.loadWeightsAsDimArray(ds_weights, "w");
fw, = sw.plotWeights(w; isBarPlot=false, label = "overall weight");
fw

# Plot generalized distances
Di = sw.loadWeightsAsDimArray(ds_weights, "Di");
figs = sw.plotDistancesPerformance(Di)
figs[1]

ds_Sij = ds_weights["Sij"];
src_names = dimnames(ds_Sij);
sources = [Array(ds_Sij[src_names[1]]), Array(ds_Sij[src_names[1]])];
Sij = DimArray(
    Array(ds_Sij),
    (Dim{Symbol(src_names[1])}(sources[1]), Dim{Symbol(src_names[2])}(sources[2])), 
    metadata = Dict(ds_Sij.attrib)
);
figs = sw.plotDistancesIndependence(Sij);
figs[1]



# TODO:
# plot distances for all combinations of diagnostics and variables
# distances_perform = sw.loadWeightsAsDimArray(ds_weights, "performance_distances");
# figs = sw.plotDistancesPerformance()




########################### 4. APPLY WEIGHTS ###########################
# simple average across model members
unweighted_means_members = sw.computeWeightedAvg(model_data_lgm.data["tas_CLIM_lgm"])

# unweighted avg across models
lgm_tas_data = sw.summarizeEnsembleMembersVector(model_data_lgm.data["tas_CLIM_lgm"], true)
unweighted_means = sw.computeWeightedAvg(lgm_tas_data)
sw.plotMeansOnMap(unweighted_means, "unweighted average LGM: tas_CLIM")
 
# weighted avg across models
weighted_means = sw.computeWeightedAvg(lgm_tas_data; weights=weights.w)
sw.plotMeansOnMap(weighted_means, "weighted means LGM: tas_CLIM")



# weighted and unweighted means should be different:
Array(weighted_means) .== Array(unweighted_means)



# apply weights

historical_model_data = sw.loadData(
    "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical",
    config_path;
    dir_per_var = false,
    subset = Dict(
        "statistics" => ["CLIM", "STD"]
    )
);



means = sw.computeWeightedAvg(data.data["tas_CLIM_lgm"])
sw.plotMeansOnMap(means, "mean LGM: tas_CLIM")




# Examples just load cmip5 or cmip6 or both
# load just cmip5 / cmip6
lgm_cmip5 = sw.loadData(
    base_path,
    config_path;
    dir_per_var = true,
    subset = Dict(
        "statistics" => ["CLIM"],
       "projects" => ["CMIP5"]
    )
);
lgm_cmip6 = sw.loadData(
    base_path,
    config_path;
    dir_per_var = true,
    subset = Dict(
        "statistics" => ["CLIM"],
       "projects" => ["CMIP6"]
    )
);
#  load both, cmip5+cmip6
lgm_cmip = sw.loadData(
    base_path,
    config_path;
    dir_per_var = true,
    subset = Dict(
        "statistics" => ["CLIM"],
        "projects" => ["CMIP"]
    )
);