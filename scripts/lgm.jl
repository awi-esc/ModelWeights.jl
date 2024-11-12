import SimilarityWeights as sw
using NCDatasets

base_path =  "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM/";
config_path = "/albedo/home/brgrus001/SimilarityWeights/configs/lgm-cmip5-cmip6";

# 1. Load model data
lgm_data = sw.loadData(
    base_path,
    config_path;
    dir_per_var = true,
    common_models_across_vars = true,
    subset = Dict(
        "statistics" => ["CLIM", "STD"],
        "data_type" => ["CMIP6"]
    )
);

# 2. Load observational data
obs_data = sw.loadData(
    "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical",
    "/albedo/home/brgrus001/SimilarityWeights/configs/recipe_configs_historical/";
    dir_per_var = true,
    isModelData = false,
    subset = Dict(
        "variables" => ["tas", "tos"], 
        "data_type" => ["ERA5"],
        #"timeranges" => ["1980-2014"]
        "aliases" => ["historical"]
    )
);

# 3. Compute weights
# if target_dir is provided within config_weights, the weights will directly 
# be saved and written to a file
config_weights = sw.ConfigWeights(
    performance = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    independence = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    ref_period = "1980-2014"#,
    # target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/"
);
weights = sw.computeWeights(lgm_data, obs_data, config_weights);

# weights can also be  saved seperately:
target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/"
target_fn = "lgm-weights.nc"
sw.saveWeights(weights, target_dir; target_fn = target_fn)


# 4. Plot weights/generalized distances
path_weights = joinpath(target_dir, "2024-11-11_12_08_lgm-weights.nc")
ds_weights = NCDataset(path_weights);

wP = sw.loadWeightsAsDimArray(ds_weights, "wP");
figs = sw.plotPerformanceWeights(wP; isBarPlot=false);
figs[1]

wI = sw.loadWeightsAsDimArray(ds_weights, "wI");
f = sw.plotWeightContributions(wI, wP)

w = sw.loadWeightsAsDimArray(ds_weights, "w");
fw = sw.plotWeights(w)

ds_Sij = ds_weights["Sij"];
src_names = dimnames(ds_Sij);
sources = [Array(ds_Sij[src_names[1]]), Array(ds_Sij[src_names[1]])];
Sij = DimArray(
    Array(ds_Sij),
    (Dim{Symbol(src_names[1])}(sources[1]), Dim{Symbol(src_names[2])}(sources[2])), 
    metadata = Dict(ds_Sij.attrib)
);
figs = sw.plotIndependenceWeights(Sij);
figs[1]

# 5. apply weights
# model weights are for ensembles, not unique ensemble members! Which means 
# that model predictions also have to be for ensembles, not members!
lgm_tas_data = sw.summarizeEnsembleMembersVector(lgm_data.data["tas_CLIM_lgm"], true)
weighted_means = sw.computeWeightedAvg(lgm_tas_data; weights=weights.w)
means = sw.computeWeightedAvg(lgm_tas_data)
sw.plotMeansOnMap(means, "unweighted average LGM: tas_CLIM")
sw.plotMeansOnMap(weighted_means, "weighted means LGM: tas_CLIM")

# weighted and unweighted means should be different:
Array(weighted_means) .== Array(means)



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
       "data_type" => ["CMIP5"]
    )
);
lgm_cmip6 = sw.loadData(
    base_path,
    config_path;
    dir_per_var = true,
    subset = Dict(
        "statistics" => ["CLIM"],
       "data_type" => ["CMIP6"]
    )
);
#  load both, cmip5+cmip6
lgm_cmip = sw.loadData(
    base_path,
    config_path;
    dir_per_var = true,
    subset = Dict(
        "statistics" => ["CLIM"],
        "data_type" => ["CMIP"]
    )
);