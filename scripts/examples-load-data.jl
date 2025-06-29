import ModelWeights as mw
import ModelWeights.Data as mwd

using NCDatasets
using DimensionalData

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

subset = Dict{String, Union{Vector{String}, mwd.Level}}(
    "statistics" => statistics, 
    "variables" => variables,
    "projects" => projects,
    "aliases" => ["lgm"], 
    "subset_shared" => mwd.MEMBER_LEVEL,
    "base_subdirs" => ["20241114"]
);
lgm_meta = mwd.loadDataFromESMValToolRecipes(
    path_data, path_recipes; subset=subset, preview=true
) 
lgm_data = mwd.loadDataFromESMValToolRecipes(
    path_data, path_recipes; subset=subset, dtype=mwd.MODEL_DATA
)

# we set subset_shared to mw.MEMBER, so model members are identical for 
# every loaded data set (variable)
model_members_lgm = Array(dims(lgm_data.models["tas_CLIM_lgm"], :member));
models_lgm = string.(unique(lgm_data.models["tos_CLIM_lgm"].properties["model_names"]))

# --------------------------------------------------------------------------- #
# 2. Model data just for historical experiment of models with lgm-experiment
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical";
path_recipes = "./configs/historical";

# 2.1 use same models (NOT on level of model members) as for the lgm 
# experiment and make sure that across variables, only shared model members 
# are loaded (subset_shared set to mw.MEMBER) since we want the exact 
# same simulations for all variables when computing weights
historical_data_lgm = mwd.loadDataFromESMValToolRecipes(
    path_data, 
    path_recipes;
    subset = Dict(
        "statistics" => statistics, 
        "variables" => variables, 
        "projects" => projects,
        "aliases" => ["historical"], 
        #"timeranges" => ["full"], 
        "models" => models_lgm
    ),
    dtype = mwd.MODEL_DATA
)
# sanity check: for all lgm models, historical experiment is loaded
models_historical = unique(historical_data_lgm.models["tos_CLIM_historical"].properties["model_names"]);
@assert models_historical == models_lgm

# function to join two DataMaps into one
data = mwd.joinDataMaps(lgm_data.models, historical_data_lgm.models)
data_members = mwd.subsetModelData(data; level = mwd.MEMBER_LEVEL)

# 2.2 Or directly load historical data of the same model MEMBERS as for lgm 
# (may be less than in 2.1, since the exact simulations now have to match with
# the lgm models, not only the models)
begin
    historical_data_lgm_members = mwd.loadDataFromESMValToolRecipes(
        path_data, 
        path_recipes;
        subset = Dict(
            "statistics" => statistics, 
            "variables" => variables, 
            "projects" => projects,
            #"aliases" => ["historical"], 
            "timeranges" => ["full"], 
            "base_subdirs" =>   ["20250211", "20250207", "20250209"],
            "models" => model_members_lgm
        ),
        dtype = mwd.MODEL_DATA
    )
end

# sanity check: all historical model members must be in lgm model members
@assert all(map(x -> x in model_members_lgm, dims(historical_data_lgm_members.models["tos_CLIM_historical"], :member)))

# sanity check: are there lgm model MEMBERS that dont appear in historical models?
members_historical = Array(dims(historical_data_lgm_members.models["tas_CLIM_historical"], :member));
filter(x -> !(x in members_historical), model_members_lgm)


# --------------- Load the model data (Version2: from yaml config file) --------------- #
# 3. Load model data for experiment lgm and historical in one run from new config file
begin
    # yaml config file already contains basic constraints for subset as defined above.
    path_config = "./configs/examples/example-lgm-historical.yml";
    subset = Dict{String, Union{Vector{String}, mwd.Level}}(
        "models" => model_members_lgm,
        "subset_shared" => mwd.MEMBER_LEVEL # applies to every loaded dataset
    );
    meta_lgm_v2 =  mwd.loadDataFromYAML(
        path_config; arg_constraint = subset, preview = true
    )
    data_lgm_v2 = mwd.loadDataFromYAML(path_config; arg_constraint = subset)
end


# -------------------- Load the observational data -------------------------- #
begin
    base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/obs/ERA5/recipe_ERA5_tas_tos_pr_20250307_174538"
    config_path = "./configs/obs"

    # aliases and timeranges don't have to match, all data will be loaded that 
    # corresponds either to aliases or to timeranges!
    obs_data = mwd.loadDataFromESMValToolRecipes(
        base_path, config_path;
        dir_per_var = false,
        subset = Dict(
            "statistics" => statistics, 
            "variables" => variables,
            #"projects" => ["ERA5"],
            "timeranges" => ["full", "1961-1990"]
            #"aliases" => ["historical", "historical2"],
        ),
        preview=false
    )
end


# --------- Load the model data Version 3: completely independent of ESMValTool recipes --------- #
base = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool" 
paths_lgm_tos = [
    joinpath(base, "LGM/recipe_cmip5_lgm_tos_20241114_150049/preproc/lgm/tos_CLIM"),
    joinpath(base, "LGM/recipe_cmip6_lgm_tos_20241114_151009/preproc/lgm/tos_CLIM")
];
paths_lgm_tas = [
    joinpath(base, "LGM/recipe_cmip6_lgm_tas_20241114_150706/preproc/lgm/tas_CLIM"),
    joinpath(base, "LGM/recipe_cmip5_lgm_tas_20241114_145900/preproc/lgm/tas_CLIM")
];
fn_format = :esmvaltool;

# returns DataMap with 1 entry
tos = mwd.loadData(paths_lgm_tos, "tos"; fn_format)
# one directory path for every dataset
paths_lgm_cmip6 = [joinpath(base, "LGM/recipe_cmip6_lgm_tos_20241114_151009/preproc/lgm/tos_CLIM"), joinpath(base, "LGM/recipe_cmip6_lgm_tas_20241114_150706/preproc/lgm/tas_CLIM")]
lgm_cmip6 = mwd.loadData(paths_lgm_cmip6, ["tos_lgm", "tas_lgm"]; fn_format)

# several directory paths for every dataset; returns DataMap with 2 entries
paths_lgm = [paths_lgm_tos, paths_lgm_tas];
lgm = mwd.loadData(paths_lgm, ["tos_lgm", "tas_lgm"]; fn_format)

shared_models = mwd.sharedModels(lgm.models, mwd.MODEL_LEVEL);
shared_members = mwd.sharedModels(lgm.models, mwd.MEMBER_LEVEL);
subset = Dict{String, Union{Vector{String}, mwd.Level}}(
    "subset_shared" => mwd.MEMBER_LEVEL
);
lgm = mwd.loadData(paths_lgm, ["tos_lgm", "tas_lgm"]; constraint=subset, fn_format)

model_members_lgm = Array(dims(lgm.models["tas_lgm"], :member));
models_lgm = string.(unique(lgm.models["tos_lgm"].properties["model_names"]))


paths_historical_tas = [
    joinpath(base, "historical/recipe_cmip5_historical_tas_20250211_094633/preproc/historical/tas_CLIM"), 
    joinpath(base, "historical/recipe_cmip6_historical_tas_20250207_080843/preproc/historical/tas_CLIM")
];
paths_historical_tos = [
    joinpath(base, "historical/recipe_cmip5_historical_tos_20250211_094633/preproc/historical/tos_CLIM"),
    joinpath(base, "historical/recipe_cmip6_historical_tos_20250209_144722/preproc/historical/tos_CLIM")
];
paths_historical = [paths_historical_tas, paths_historical_tos];
subset = Dict{String, Union{Vector{String}, mwd.Level}}(
    "models" => model_members_lgm,
    "subset_shared" => mwd.MEMBER_LEVEL # applies to every loaded dataset
);
historical = mwd.loadData(
    paths_historical, ["tas_historical", "tos_historical"]; fn_format, constraint=subset
)

# to load data as YAXArrays from files directly (not from all files within directories), use loadPreprocData
paths = vcat(mwd.collectNCFilePaths.(paths_lgm_tas)...)
data = mwd.loadPreprocData(
    paths, fn_format; dtype=mwd.MODEL_DATA, meta_info=Dict("variable" => "tas")
)
# same data but a ClimateData instance is loaded 
data = mwd.loadData(paths, "tas"; fn_format, meta_data = Dict("variable" => "tas"))