import ModelWeights as mw
using NCDatasets
using DimensionalData

########################### 1. LOADING DATA ###########################
# Model data just for lgm-experiment from ESMValTool recipes
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/lgm-cmip5-cmip6";
lgm_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes;
    dir_per_var = true, # default: true
    is_model_data = true, # default: true 
    only_shared_models = true, # default: false
    subset = mw.Constraint(
        statistics = ["CLIM"],
        variables = ["tas", "tos"],
        aliases = ["lgm"],
        projects = ["CMIP5", "CMIP6"],
        models = Vector{String}(),
        subdirs = ["20241114"]
    ),
    preview = false # default value is false
);
# we set only_shared_models to true, so model members are identical for every 
# loaded data set 
model_members_lgm = Array(dims(lgm_data[1].data, :member))
models_lgm = unique(lgm_data[1].data.metadata["model_names"])

# Model data for historical experiment of models with lgm-experiment from above
base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/";
config_path = "/albedo/home/brgrus001/ModelWeights/configs/historical";
historical_data = mw.loadDataFromESMValToolConfigs(
    base_path, config_path;
    only_shared_models = true,
    subset = mw.Constraint(
        statistics = ["CLIM"],
        variables = ["tas", "tos"],
        aliases = ["historical"],
        timeranges = ["full"],
        models = model_members_lgm,
        #models = models_lgm,
        # if dir_per_var=true names of data subdirs must contain any of:
        subdirs = ["20241121", "20241118"]
    )
);

# Load model data for experiment lgm and historical in one run from new config file
path_config = "/albedo/home/brgrus001/ModelWeights/configs/examples/example-lgm-historical.yml";
model_data = mw.loadDataFromYAML(
    path_config;
    dir_per_var = true, # true is default value
    is_model_data = true, # true is default value
    only_shared_models = true,
    # subset = mw.Constraint(
    #     projects = ["CMIP5"]
    #     "models" => model_members_lgm
    # ),
    preview = false
);


# Load the observational data
base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_obs_historical_20241119_124434";
config_path = "/albedo/home/brgrus001/ModelWeights/configs/historical_obs";
obs_data = mw.loadDataFromESMValToolConfigs(
    base_path, config_path;
    dir_per_var = false,
    is_model_data = false,
    subset = mw.Constraint(
        statistics = ["CLIM"],
        variables = ["tas", "tos"],
        aliases = ["historical", "historical2"],
        projects = ["ERA5"], # for observational data default value is ["ERA5"]
        timeranges = ["1980-2014"]
    )
);


########################### 2. COMPUTATION WEIGHTS ###########################
# if target_dir is provided within config_weights, the weights will directly 
# be saved and written to a file
config_weights = mw.ConfigWeights(
    performance = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    independence = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    # ref_period can refer to either an alias or a timerange:
    ref_period = "historical", # alias
    #ref_period = "full", # timerange
    target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/"
);
weights = mw.computeWeights(historical_data, obs_data, config_weights);
#weights = mw.computeWeights(model_data, obs_data, config_weights);

# weights can also be  saved separately:
target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/";
target_fn = "weights-lgm-models.nc";
mw.saveWeightsAsNCFile(weights, target_dir; target_fn = target_fn)
target_path = joinpath(target_dir, "weights-lgm-models.jld2");
mw.saveWeightsAsJuliaObj(weights, target_path)


#  Load weights from/as Julia object
weights = mw.loadWeightsFromJLD2(target_path)


########################### 3. PLOTTING ###########################
# Plot weights/generalized distances
path_weights = joinpath(target_dir, target_fn);
ds_weights = NCDataset(path_weights);


ds_weights = NCDataset("/albedo/work/projects/p_pool_clim_data/britta/weights/weights_2024-11-27_09_22.nc")
# Plot performance weights
wP = mw.loadWeightsAsDimArray(ds_weights, "wP");
fig_wP, = mw.plotWeights(wP; isBarPlot=false, label="performance weight");
fig_wP

# Plot independence weights
wI = mw.loadWeightsAsDimArray(ds_weights, "wI");
fig_wI, = mw.plotWeights(wI; isBarPlot=false, label="independence weight");
fig_wI

# Plot performance and independence weights together
f = mw.plotWeightContributions(wI, wP)

# Plot overall weights
w = mw.loadWeightsAsDimArray(ds_weights, "w");
fw, = mw.plotWeights(w; isBarPlot=false, label = "overall weight");
fw

# Plot generalized distances
Di = mw.loadWeightsAsDimArray(ds_weights, "Di");
figs = mw.plotDistancesPerformance(Di)
figs[1]

ds_Sij = ds_weights["Sij"];
src_names = dimnames(ds_Sij);
sources = [Array(ds_Sij[src_names[1]]), Array(ds_Sij[src_names[1]])];
Sij = DimArray(
    Array(ds_Sij),
    (Dim{Symbol(src_names[1])}(sources[1]), Dim{Symbol(src_names[2])}(sources[2])), 
    metadata = Dict(ds_Sij.attrib)
);
figs = mw.plotDistancesIndependence(Sij);
figs[1]



# TODO:
# plot distances for all combinations of diagnostics and variables
# distances_perform = mw.loadWeightsAsDimArray(ds_weights, "performance_distances");
# figs = mw.plotDistancesPerformance()




########################### 4. APPLY WEIGHTS ###########################
# simple average across model members
unweighted_means_members = mw.computeWeightedAvg(model_data_lgm.data["tas_CLIM_lgm"])

# unweighted avg across models
lgm_tas_data = mw.summarizeEnsembleMembersVector(model_data_lgm.data["tas_CLIM_lgm"], true)
unweighted_means = mw.computeWeightedAvg(lgm_tas_data)
mw.plotMeansOnMap(unweighted_means, "unweighted average LGM: tas_CLIM")
 
# weighted avg across models
weighted_means = mw.computeWeightedAvg(lgm_tas_data; weights=weights.w)
mw.plotMeansOnMap(weighted_means, "weighted means LGM: tas_CLIM")



# weighted and unweighted means should be different:
Array(weighted_means) .== Array(unweighted_means)



# apply weights

historical_model_data = mw.loadData(
    "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical",
    config_path;
    dir_per_var = false,
    subset = Dict(
        "statistics" => ["CLIM", "STD"]
    )
);



means = mw.computeWeightedAvg(data.data["tas_CLIM_lgm"])
mw.plotMeansOnMap(means, "mean LGM: tas_CLIM")




# Examples just load cmip5 or cmip6 or both
# load just cmip5 / cmip6
lgm_cmip5 = mw.loadData(
    base_path,
    config_path;
    dir_per_var = true,
    subset = Dict(
        "statistics" => ["CLIM"],
       "projects" => ["CMIP5"]
    )
);
lgm_cmip6 = mw.loadData(
    base_path,
    config_path;
    dir_per_var = true,
    subset = Dict(
        "statistics" => ["CLIM"],
       "projects" => ["CMIP6"]
    )
);
#  load both, cmip5+cmip6
lgm_cmip = mw.loadData(
    base_path,
    config_path;
    dir_per_var = true,
    subset = Dict(
        "statistics" => ["CLIM"],
        "projects" => ["CMIP"]
    )
);