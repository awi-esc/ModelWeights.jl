import ModelWeights as mw
using NCDatasets
using DimensionalData
using Setfield

# --------------------------- Set configurations --------------------------- #
begin
    dir_per_var = true;
    is_model_data = true;
    statistics = ["CLIM"];
    variables = ["tas", "tos"];
    projects = ["CMIP5", "CMIP6"];
end
# --------- Load the model data Version 1: from ESMValTool recipes --------- #
# 1. Model data just for lgm-experiment from ESMValTool recipes
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/lgm-cmip5-cmip6";

subset = Dict{String, Union{Vector{String}, mw.LEVEL}}(
    "statistics" => statistics, 
    "variables" => variables, 
    "projects" => projects,
    "aliases" => ["lgm"], 
    "subdirs" => ["20241114"],
    "level_shared_models" => mw.MEMBER
);
lgm_meta = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes; dir_per_var, is_model_data, subset = subset, 
    preview = true
);
lgm_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes; dir_per_var, is_model_data, subset = subset,
    preview = false
);

# we set level_shared_models to mw.MEMBER, so model members are identical for 
# every loaded data set (variable)
model_members_lgm = Array(dims(lgm_data["tas_CLIM_lgm"].data, :member));
models_lgm = unique(lgm_data["tos_CLIM_lgm"].data.metadata["model_names"]);

# --------------------------------------------------------------------------- #
# 2. Model data just for historical experiment of models with lgm-experiment
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/historical";

# 2.1 use same models (NOT on level of model members) as for the lgm 
# experiment and make sure that across variables, only shared model members 
# are loaded (level_shared_models set to mw.MEMBER) since we want the exact 
# same simulations for all variables when computing weights
historical_data_lgm_models = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes;
    subset = Dict(
        "statistics" => statistics, 
        "variables" => variables, 
        "projects" => projects,
        "aliases" => ["historical"], 
        "timeranges" => ["full"], 
        "subdirs" =>  ["20241121", "20241118"],
        "models" => models_lgm
    ),
    preview = false
);
# sanity check: for all lgm models, historical experiment is loaded
models_historical = unique(historical_data_lgm_models["tos_CLIM_historical"].data.metadata["model_names"]);
@assert models_historical == models_lgm


# 2.2 Or load historical data of the same model MEMBERS as for lgm 
# (may be less than in 2.1, since the exact simulations now have to match with
# the lgm models, not only the models)
begin
    historical_data_lgm_members = mw.loadDataFromESMValToolConfigs(
        path_data, path_recipes;
        subset = Dict(
            "statistics" => statistics, 
            "variables" => variables, 
            "projects" => projects,
            "aliases" => ["historical"], 
            "timeranges" => ["full"], 
            "subdirs" =>  ["20241121", "20241118"],
            "models" => model_members_lgm
        ),
        preview = false
    );
end

# sanity check: all historical model members must be in lgm model members
all(map(x -> x in model_members_lgm, dims(historical_data_lgm_members["tos_CLIM_historical"].data, :member)))

# sanity check: are there lgm model members that dont appear in historical models?
members_historical = Array(dims(historical_data_lgm_members["tas_CLIM_historical"].data, :member));
filter(x -> !(x in members_historical), model_members_lgm)


# --------------- Load the model data (Version2: from yaml config file) --------------- #
# 3. Load model data for experiment lgm and historical in one run from new config file
begin
    # yaml config file already contains basic constraints for subset as defined above.
    path_config = "/albedo/home/brgrus001/ModelWeights/configs/examples/example-lgm-historical.yml";
    subset = Dict{String, Union{Vector{String}, mw.LEVEL}}(
        "statistics" => statistics, 
        "variables" => variables, 
        "projects" => projects,
        "aliases" => ["historical"],
        "models" => models_lgm
    );
    hist_data_lgm_v2_meta =  mw.loadDataFromYAML(
        path_config; 
        is_model_data,
        subset = subset,
        preview = true
    );
    hist_data_lgm_v2 = mw.loadDataFromYAML(
        path_config; 
        is_model_data,
        subset = subset,
        preview = false
    );
end

# TODO: make some checks that same data is loaded with both versions
# begin
#     missingOrEqual(arr1, arr2) = all(x -> ismissing(x) || x==1, Array(arr1) .== Array(arr2));
#     for v in ["tas", "tos"]
#         data1 = historical_data_lgm_members[v * "_CLIM_historical"].data;
#         data2 = historical_lgm_data_config[v * "_CLIM_historical"].data;
#         @assert missingOrEqual(data1, data2)
#     end
# end



# -------------------- Load the observational data -------------------------- #
begin
    base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_obs_historical_tas_tos";
    config_path = "/albedo/home/brgrus001/ModelWeights/configs/historical_obs";
    # aliases and timeranges don't have to match, all data will be loaded that 
    # corresponds either to aliases or to timeranges!
    obs_data = mw.loadDataFromESMValToolConfigs(
        base_path, config_path;
        dir_per_var = false,
        is_model_data = false,
        subset = Dict(
            "statistics" => statistics, 
            "variables" => variables,
            "aliases" => ["historical", "historical2"],
            "projects" => ["ERA5"], # for observational data default value is ["ERA5"]
            "timeranges" => ["1980-2014"]
        )
    );
end