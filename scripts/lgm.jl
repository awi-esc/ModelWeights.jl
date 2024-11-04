import SimilarityWeights as sw

path_lgm_data =  "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM/";

#target_dir = "recipe_cmip6_lgm_tas" #"recipe_lgm_20241021_061102";
config_path = "/albedo/home/brgrus001/SimilarityWeights/configs/lgm-cmip5-cmip6";

base_path = path_lgm_data;# * target_dir;

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

config_weights = sw.ConfigWeights(
    performance = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    independence = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    ref_period = "1980-2014"
);


weights = sw.getOverallWeights(lgm_data, obs_data, config_weights)
Array(weights.wP)


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