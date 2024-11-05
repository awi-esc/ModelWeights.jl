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
        "statistics" => ["CLIM"],
        "data_type" => ["CMIP"]
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
        "aliases" => ["historical"]
    )
)

# 3. Compute weights
config_weights = sw.ConfigWeights(
    performance = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    independence = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    ref_period = "1980-2014"#,
    #target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/"
    );
    
# if target_dir is provided within config_weights, the weights will directly 
# be saved and written to a file
weights = sw.getOverallWeights(lgm_data, obs_data, config_weights);

# or save weights seperately afterwards
target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/"
target_fn = "lgm-weights.nc"
sw.saveWeights(weights, target_dir; target_fn = target_fn)


# 4. Plot weights
path_weights = joinpath(target_dir, "lgm-weights.nc")
weights = NCDataset(path_weights)

wP = sw.loadWeightsAsDimArray(weights, "wP")
figs = sw.plotPerformanceWeights(wP; isBarPlot=false)

wI = sw.loadWeightsAsDimArray(weights, "wI")
f = sw.plotWeightContributions(wI, wP)


w = sw.loadWeightsAsDimArray(weights, "w")
fw = sw.plotWeights(w)
#TODO: also save generalized distances besides normalized weights
#Sij = sw.loadWeightsAsDimArray(weights, "Sij")
#figs = sw.plotIndependenceWeights(wI)







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