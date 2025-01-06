import ModelWeights as mw
using NCDatasets
using DimensionalData
using Setfield

########################### 1. LOADING DATA ###########################

# ------------------------ Load the model data ------------------------------ #
# Model data just for lgm-experiment from ESMValTool recipes
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/lgm-cmip5-cmip6";

dir_per_var = true;
is_model_data = true;
only_shared_models = true;
statistics = ["CLIM"];
variables = ["tas", "tos"];
projects = ["CMIP5", "CMIP6"];
subset = mw.Constraint(
    statistics = statistics, variables = variables, projects = projects,
    aliases = ["lgm"], subdirs = ["20241114"]
);

lgm_meta = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes; dir_per_var, is_model_data, only_shared_models,
    subset, preview = true
);
lgm_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes; dir_per_var, is_model_data, only_shared_models,
    subset, preview = false
);

# we set only_shared_models to true, so model members are identical for every 
# loaded data set 
model_members_lgm = Array(dims(lgm_data["tos_CLIM_lgm"].data, :member));
models_lgm = unique(lgm_data["tos_CLIM_lgm"].data.metadata["model_names"]);

# --------------------------------------------------------------------------- #
# Model data for historical experiment of models with lgm-experiment from above
# Version1: load data from ESMValToolConfigs and subset to lgm_models
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/historical";

# use same big models (NOT on level of model members)
model_data0 = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes;
    only_shared_models = true,
    subset = mw.Constraint(
        statistics = statistics, variables = variables, projects = projects,
        aliases = ["historical"], 
        timeranges=["full"], 
        subdirs =  ["20241121", "20241118"],
        models = models_lgm,
    ),
    preview = false
);
# sanity check: all lgm models are in historical data?
models_historical = unique(model_data0["tos_CLIM_historical"].data.metadata["model_names"]);
@assert models_historical == models_lgm


# load historical data of the same model members as in lgm data
begin
    model_historical = mw.loadDataFromESMValToolConfigs(
        path_data, path_recipes;
        only_shared_models = true,
        subset = mw.Constraint(
            statistics = statistics, variables = variables, projects = projects,
            aliases = ["historical"], 
            timeranges=["full"], 
            subdirs =  ["20241121", "20241118"],
            models = model_members_lgm,
        ),
        preview = false
    );
end

# sanity check: all lgm models are in historical data?
model_members_historical = Array(dims(model_historical["tas_CLIM_historical"].data, :member));
filter(x -> !(x in model_members_historical), model_members_lgm)

# Version2: Load model data for experiment lgm and historical in one run from new config file
begin
    path_config = "/albedo/home/brgrus001/ModelWeights/configs/examples/example-lgm-historical.yml";
    # yaml config file already contains basic constraints for subset as defined above.
    model_historical_lgm = mw.loadDataFromYAML(path_config; dir_per_var, is_model_data,
        only_shared_models = true,
        subset = mw.Constraint(models = model_members_lgm),
        # subset = mw.Constraint(models = models_lgm),
        preview = false
    );
end

begin
    # some checks that same data is loaded with both versions
    missingOrEqual(arr1, arr2) = all(x -> ismissing(x) || x==1, Array(arr1) .== Array(arr2));
    for v in ["tas", "tos"]
        data1 = model_historical[v * "_CLIM_historical"].data;
        data2 = model_historical_lgm[v * "_CLIM_historical"].data;
        @assert missingOrEqual(data1, data2)
    end
end

# -------------------- Load the observational data -------------------------- #
begin
    base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_obs_historical_20241119_124434";
    config_path = "/albedo/home/brgrus001/ModelWeights/configs/historical_obs";
    # aliases and timeranges don't have to match, all data will be loaded that 
    # corresponds either to aliases or to timeranges!
    obs_data = mw.loadDataFromESMValToolConfigs(
        base_path, config_path;
        dir_per_var = false,
        is_model_data = false,
        subset = mw.Constraint(statistics = statistics, variables = variables,
            aliases = ["historical", "historical2"],
            projects = ["ERA5"], # for observational data default value is ["ERA5"]
            timeranges = ["1980-2014"]
        )
    );
end

########################### 2. COMPUTATION WEIGHTS ###########################
# if target_dir is provided within config_weights, the weights will directly 
# be saved and written to a file
begin
    target_dir = "/albedo/work/projects/p_pool_clim_data/britta/weights/";
    target_fn = "weights-lgm.nc";
    # Note: ref_perform_weights can refer to either an alias or a timerange!
    config_weights = mw.ConfigWeights(
        performance = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
        independence = Dict("tas_CLIM"=>1, "tos_CLIM"=>1),
        sigma_independence = 0.5,
        sigma_performance = 0.5,
        ref_perform_weights = "historical", # alias
        #ref_perform_weights = "full", # timerange
        ref_indep_weights = "historical", # alias
        target_path = joinpath(target_dir, target_fn)
    );
    weights = mw.computeWeights(collect(values(model_historical)), collect(values(obs_data)), config_weights);
end

begin
    # Weights for lgm experiment, with performance weights wrt historical data, but independence 
    # weights wrt lgm experiments
    config_weights_lgm = @set config_weights.ref_indep_weights = "lgm";
    weights_lgm = mw.computeWeights(collect(values(model_historical_lgm)), collect(values(obs_data)), config_weights_lgm);
end

# weights can also be  saved separately (as julia obj or .nc file):
target_nc = mw.saveWeightsAsNCFile(weights, joinpath(target_dir, "weights-lgm.nc"));
target_jld2 = mw.saveWeightsAsJuliaObj(weights, joinpath(target_dir, "weights-lgm.jld2"));
#  Load weights as Julia object / NCDataset
weights_from_disk = mw.loadWeightsFromJLD2(target_jld2);
ds_weights = NCDataset(target_nc);

@assert target_jld2 == weights_from_disk.config.target_path

########################### 3. PLOTTING ###########################
# weights = weights_from_disk
# if weights loaded from .nc file instead, you can use the following function to load
# respective weights:
# wP = mw.loadWeightsAsDimArray(ds_weights, "wP");
# Plot weights/generalized distances

# Plot performance weights
fig_wP, = mw.plotWeights(weights.wP; is_bar_plot = false, label = "performance weight");
fig_wP

# Plot independence weights
fig_wI, = mw.plotWeights(weights.wI; is_bar_plot=false, label="independence weight");
fig_wI

# Plot performance and independence weights together
f = mw.plotWeightContributions(wI, wP)

# Plot overall weights
fw, = mw.plotWeights(weights.w; is_bar_plot=false, label = "overall weight");
fw

# Plot generalized distances
figs = mw.plotDistancesPerformance(weights.Di)
figs[1]

ds_Sij = ds_weights["Sij"];
src_names = dimnames(ds_Sij);
sources = [Array(ds_Sij[src_names[1]]), Array(ds_Sij[src_names[1]])];
Sij = DimArray(
    Array(ds_Sij),
    (Dim{Symbol(src_names[1])}(sources[1]), Dim{Symbol(src_names[2])}(sources[2])), 
    metadata = Dict(ds_Sij.attrib)
);
# Note for weights.Sij dimension names are 'model1' and 'model2', but fn expects twice 'model'
figs = mw.plotDistancesIndependence(Sij);
figs[1]


figs_performances = mw.plotDistancesPerformance(weights.performance_distances);
figs_performances[1]
figs_performances[2]

# TODO: fix plotting for independence distances, dimension names are 'model1' and 'model2', but fn expects twice 'model'
#figs = mw.plotDistancesIndependence(weights.independence_distances);


########################### 4. APPLY WEIGHTS ###########################
tas_data = model_historical["tas_CLIM_historical"].data;
# The function applyWeights computes the weighted avg based on the average 
# across the respective members of each model, if the data is given on the
# level of model members. It further checks that weights were computed for the 
# same models as in the provided data, and possibly takes a subset of the data 
# or models (in that case weights are renormalized)
weighted_avg_models = mw.applyWeights(tas_data, weights.w);
# if you want to compute the weighted avg based on the weights for each model 
# on the level of model members,use the function computeWeightedAvg directly:
weighted_avg_members = mw.computeWeightedAvg(tas_data; weights=weights.w_members);


unweighted_avg = mw.computeWeightedAvg(tas_data; use_members_equal_weights = false);
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
mw.plotMeansOnMap(unweighted_means, "unweighted average LGM: tas_CLIM")
 
# weighted avg across models
weighted_means = mw.computeWeightedAvg(lgm_tas_data; weights=weights.w)
mw.plotMeansOnMap(weighted_means, "weighted means LGM: tas_CLIM")

# weighted and unweighted means should be different:
Array(weighted_means) .== Array(unweighted_means)



