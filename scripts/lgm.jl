import SimilarityWeights as sw
using NCDatasets
using DimensionalData


# 1. Load model data
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
        
model_members_lgm = Array(dims(first(values(model_data_lgm.data)), :member))
models_lgm  = unique(first(values(model_data_lgm.data)).metadata["model_names"])


base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/";
config_path = "/albedo/home/brgrus001/SimilarityWeights/configs/historical";
model_data_historical = sw.loadData(
    base_path, config_path;
    dir_per_var = true,
    common_models_across_vars = true,
    subset = Dict(
        "statistic" => ["CLIM"],
        "variable" => ["tas", "tos"],
        "alias" => ["historical"],
        "projects" => ["CMIP5", "CMIP6"],
        "models" => model_members_lgm,
        "subdirs" => ["20241121", "20241118"] # if dir_per_var is true only subdirs containing any are considered
    )
);
 
#model_data_historical = sw.getCommonModelsAcrossVars(model_data_historical, :member)


# 2. Load observational data
obs_data = sw.loadData(
    "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_obs_historical_20241119_124434",
    "/albedo/home/brgrus001/SimilarityWeights/configs/historical_obs";
    dir_per_var = false,
    isModelData = false,
    subset = Dict(
        "statistic" => ["CLIM"],
        "variable" => ["tas", "tos"],
        "alias" => ["historical"],
        "projects" => ["ERA5"] # default value (same as not setting it)
        #"timeranges" => ["1980-2014"]
    )
);

# 3. Compute weights
# if target_dir is provided within config_weights, the weights will directly 
# be saved and written to a file
config_weights = sw.ConfigWeights(
    performance = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    independence = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    sigma_independence = 0.5,
    sigma_performance = 0.5#,
    #ref_period = "1980-2014"#,
    # target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/"
);
weights = sw.computeWeights(model_data_historical, obs_data, config_weights);

# weights can also be  saved seperately:
target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/"
target_fn = "weights-lgm-models.nc"
sw.saveWeights(weights, target_dir; target_fn = target_fn)


# 4. Plot weights/generalized distances
path_weights = joinpath(target_dir, target_fn)
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