import ModelWeights as mw
import ModelWeights.Data as mwd

using DimensionalData
using NCDatasets
using YAXArrays

# --------------------------- Set general configurations --------------------------- #
filename_format = :esmvaltool;
dtype = "cmip"

# --------- Load the model data Version 1: from ESMValTool recipes --------- #
dir_per_var = true;
statistics = ["CLIM"];
variables = ["tas", "tos"];
mips = ["CMIP5", "CMIP6"];

# 1. Model data just for lgm-experiment from ESMValTool recipes
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/lgm-cmip5-cmip6";

constraint = Dict{String, Union{Vector{String}, String}}(
    "statistics" => statistics, 
    "variables" => variables,
    "mips" => mips,
    "aliases" => ["lgm"], 
    "level_shared" => "member",
    "base_subdirs" => ["20241114"]
);
lgm_meta = mw.defineDataMap(
    path_data, path_recipes, :esmvaltool_recipes; constraint, preview = true
) 
lgm_data = mw.defineDataMap(
    path_data, path_recipes, :esmvaltool_recipes; constraint, dtype = "cmip"
)

# we set level_shared to mw.MEMBER, so model members are identical for 
# every loaded data set (variable)
model_members_lgm = Array(dims(lgm_data["tas_CLIM_lgm"], :member));
models_lgm = mwd.modelsFromMemberIDs(model_members_lgm; uniq = true)

# --------------------------------------------------------------------------- #
# 2. Model data just for historical experiment of models with lgm-experiment

# 2.1 use same models (NOT on level of model members) as for the lgm 
# experiment and make sure that across variables, only shared model members 
# are loaded (level_shared set to mw.MEMBER) since we want the exact 
# same simulations for all variables when computing weights
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical";
path_recipes = "./configs/historical";
historical_data_lgm = mw.defineDataMap(
    path_data, 
    path_recipes,
    :esmvaltool_recipes;
    constraint = Dict(
        "statistics" => statistics, 
        "variables" => variables, 
        "mips" => mips,
        "aliases" => ["historical"],
        "models" => models_lgm
    ),
    dtype = "cmip"
)

# sanity check: for all lgm models, historical experiment is loaded
models_historical = mwd.modelsFromMemberIDs(
    val(historical_data_lgm["tos_CLIM_historical"].member); uniq = true
)
@assert models_historical == models_lgm

# function to join two DataMaps into one
data = mwd.joinDataMaps(lgm_data, historical_data_lgm)
data_members = mwd.subsetModelData(data, "member")

# 2.2 Or directly load historical data of the same model MEMBERS as for lgm 
# (may be less than in 2.1, since the exact simulations now have to match with
# the lgm models, not only the models)
begin
    historical_data_lgm_members = mwd.defineDataMap(
        path_data, 
        path_recipes,
        :esmvaltool_recipes;
        constraint = Dict(
            "statistics" => statistics, 
            "variables" => variables, 
            "mips" => mips,
            "timeranges" => ["full"], 
            "base_subdirs" =>   ["20250211", "20250207", "20250209"],
            "models" => model_members_lgm
        ),
        dtype = "cmip"
    )
end

# sanity check: all historical model members must be in lgm model members
@assert all(map(x -> x in model_members_lgm, dims(historical_data_lgm_members["tos_CLIM_historical"], :member)))

# sanity check: are there lgm model MEMBERS that dont appear in historical models?
members_historical = Array(dims(historical_data_lgm_members["tas_CLIM_historical"], :member));
filter(x -> !(x in members_historical), model_members_lgm)


# --------------- Load the model data (Version2: from yaml config file) --------------- #
# 3. Load model data for experiment lgm and historical in one run from new config file
begin
    # yaml config file already contains basic constraints for subset as defined above.
    path_config = "./configs/examples/example-lgm-historical.yml";
    constraint = Dict{String, Union{Vector{String}, Symbol}}(
        "models" => model_members_lgm,
        "level_shared" => :member # applies to every variable
    );
    meta_lgm_v2 =  mw.defineDataMap(path_config; constraint, preview=true)
    data_lgm_v2 = mw.defineDataMap(path_config; constraint)
end


# -------------------- Load the observational data -------------------------- #
begin
    base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/obs/ERA5/recipe_ERA5_tas_tos_pr_20250307_174538"
    config_path = "./configs/obs"

    # aliases and timeranges don't have to match, all data will be loaded that 
    # corresponds either to aliases or to timeranges!
    obs_data = mw.defineDataMap(
        base_path, 
        config_path,
        :esmvaltool_recipes;
        dir_per_var = false,
        constraint = Dict(
            "statistics" => statistics, 
            "variables" => variables,
            "timeranges" => ["full", "1961-1990"]
        ),
        preview = false
    )
end


# --------- Load the model data Version 3: completely independent of ESMValTool recipes --------- #
# 1. most basic case, where some available data (YAXArrays) are loaded into a DataMap
longitudes =  [12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5]
latitudes = [-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5]

v1 = YAXArray((Dim{:lat}(latitudes), Dim{:lon}(longitudes)), rand(7, 9))
v2 = YAXArray((Dim{:lat}(latitudes), Dim{:lon}(longitudes)), rand(7, 9)) 

variables = mw.defineDataMap([v1, v2], ["ESM1", "ESM2"])

# 2. Load from data directories (defineDataMap(Vector{Vector{String}}, String))
base = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool" 

# 2.1 just CMIP5, single dataset
paths_lgm_tos = [
    joinpath(base, "LGM/recipe_cmip5_lgm_tos_20241114_150049/preproc/lgm/tos_CLIM")
];
preview_tos = mwd.previewDataMap(paths_lgm_tos, "tos"; filename_format = :esmvaltool_cmip5)
tos = mwd.defineDataMap(paths_lgm_tos, "tos"; filename_format = :esmvaltool_cmip5)

# 2.2 just CMIP6, single dataset
paths_lgm_tas = [
    joinpath(base, "LGM/recipe_cmip6_lgm_tas_20241114_150706/preproc/lgm/tas_CLIM")
];
preview_tas = mwd.previewDataMap(paths_lgm_tas, "tas"; filename_format = :esmvaltool_cmip6)
tas = mwd.defineDataMap(paths_lgm_tas, "tas"; filename_format = :esmvaltool_cmip6)
# also works with filename_format = :esmvaltool:
tas = mwd.defineDataMap(paths_lgm_tas, "tas"; filename_format = :esmvaltool)

# 2.3 single dataset, but for each CMIP5 and CMIP6 data is loaded
paths_lgm_tos = [
    joinpath(base, "LGM/recipe_cmip5_lgm_tos_20241114_150049/preproc/lgm/tos_CLIM"),
    joinpath(base, "LGM/recipe_cmip6_lgm_tos_20241114_151009/preproc/lgm/tos_CLIM")
];
paths_lgm_tas = [
    joinpath(base, "LGM/recipe_cmip6_lgm_tas_20241114_150706/preproc/lgm/tas_CLIM"),
    joinpath(base, "LGM/recipe_cmip5_lgm_tas_20241114_145900/preproc/lgm/tas_CLIM")
];
# each returns a DataMap with 1 entry
preview_tas = mwd.previewDataMap(paths_lgm_tas, "tas"; filename_format=:esmvaltool)
tas = mwd.defineDataMap(paths_lgm_tas, "tas"; filename_format=:esmvaltool)
preview_tos = mwd.previewDataMap(paths_lgm_tos, "tos"; filename_format=:esmvaltool)
tos = mwd.defineDataMap(paths_lgm_tos, "tos"; filename_format=:esmvaltool) 


# 2.4 one directory path for every dataset (defineDataMap(Vector{String}, Vector{String}))
paths_lgm_cmip6 = [joinpath(base, "LGM/recipe_cmip6_lgm_tos_20241114_151009/preproc/lgm/tos_CLIM"), joinpath(base, "LGM/recipe_cmip6_lgm_tas_20241114_150706/preproc/lgm/tas_CLIM")]
preview_lgm_cmip6 = mwd.previewDataMap(paths_lgm_cmip6, ["tos_lgm", "tas_lgm"]; filename_format=:esmvaltool_cmip6)
lgm_cmip6 = mw.defineDataMap(paths_lgm_cmip6, ["tos_lgm", "tas_lgm"]; filename_format=:esmvaltool_cmip6)
# also works for (less specific) filename_format=:esmvaltool
lgm_cmip6 = mw.defineDataMap(paths_lgm_cmip6, ["tos_lgm", "tas_lgm"]; filename_format=:esmvaltool)

# same with cmip5
paths_lgm_cmip5 = [
    joinpath(base,"LGM/recipe_cmip5_lgm_tos_20241114_150049/preproc/lgm/tos_CLIM"), 
    joinpath(base, "LGM/recipe_cmip5_lgm_tas_20241114_145900/preproc/lgm/tas_CLIM")
]
preview_lgm_cmip5 = mwd.previewDataMap(paths_lgm_cmip5, ["tos_lgm", "tas_lgm"]; filename_format=:esmvaltool_cmip5)
lgm_cmip5 = mw.defineDataMap(paths_lgm_cmip5, ["tos_lgm", "tas_lgm"]; filename_format=:esmvaltool_cmip5)
lgm_cmip5 = mw.defineDataMap(paths_lgm_cmip5, ["tos_lgm", "tas_lgm"]; filename_format=:esmvaltool)

# several directory paths for every dataset; returns DataMap with 2 entries (defineDataMap(Vector{Vector{String}}, Vector{String}))
paths_lgm = [paths_lgm_tas, paths_lgm_tos];
preview_lgm = mwd.previewDataMap(paths_lgm, ["tas", "tos"]; filename_format=:esmvaltool)
lgm_cmip = mwd.defineDataMap(paths_lgm, ["tas", "tos"]; filename_format=:esmvaltool)
# add constraints
preview_lgm_cmip5 = mwd.previewDataMap(
    paths_lgm, ["tas", "tos"]; filename_format=:esmvaltool, constraint = Dict(:mips => ["CMIP5"])
)
lgm_cmip5 = mw.defineDataMap(
    paths_lgm, ["tas", "tos"]; filename_format=:esmvaltool, constraint = Dict(:mips => ["CMIP5"])
) 
preview_lgm_cmip6 = mwd.previewDataMap(
    paths_lgm, ["tas", "tos"]; filename_format=:esmvaltool, constraint = Dict(:mips => ["CMIP6"])
)
lgm_cmip6 = mw.defineDataMap(
    paths_lgm, ["tas", "tos"]; filename_format=:esmvaltool, constraint = Dict(:mips => ["CMIP6"])
)

shared_models = mwd.sharedModels(lgm, :model)
shared_members = mwd.sharedModels(lgm, :member)
model_members_lgm = Array(dims(lgm["tas"], :member))
model_names =  mwd.modelsFromMemberIDs(model_members_lgm; uniq = true)
model_names =  mwd.modelsFromMemberIDs(model_members_lgm; uniq = false)

# Adding constraints
constraint = Dict(:variables => ["tas"])
lgm_cmip6_tas = mw.defineDataMap(
    paths_lgm_cmip6, ["tos_lgm", "tas_lgm"]; filename_format=:esmvaltool_cmip6, constraint
)





# Historical Data
paths_historical_tas = [
    joinpath(base, "historical/recipe_cmip5_historical_tas_20250211_094633/preproc/historical/tas_CLIM"), 
    joinpath(base, "historical/recipe_cmip6_historical_tas_20250207_080843/preproc/historical/tas_CLIM")
];
paths_historical_tos = [
    joinpath(base, "historical/recipe_cmip5_historical_tos_20250211_094633/preproc/historical/tos_CLIM"),
    joinpath(base, "historical/recipe_cmip6_historical_tos_20250209_144722/preproc/historical/tos_CLIM")
];

# TODO: fix constraints -> if just a single dataset, level should not have an influence!
paths_historical = [paths_historical_tas, paths_historical_tos];
constraint = Dict(:models => model_members_lgm)
preview_historical = mwd.previewDataMap(
    paths_historical_tas, 
    "tas"; 
    #level = :member, 
    filename_format = :esmvaltool,
    constraint
)
historical = mwd.defineDataMap(
    paths_historical_tos, 
    "tos"; 
    level = :member, 
    filename_format = :esmvaltool, 
    #constraint
)




# to load data as YAXArrays from files directly (not from all files within directories), use loadPreprocData
paths_tas = vcat(mwd.collectNCFilePaths.(paths_lgm_tas)...)
paths_tos = vcat(mwd.collectNCFilePaths.(paths_lgm_tos)...)
data = mwd.loadPreprocData(paths_tas, filename_format; dtype="cmip")
# when cmip is not defined, default names are used for models
data = mwd.loadPreprocData(paths_tas, filename_format)
# same data but a DataMap instance is returned
data = mwd.loadDataMapCore([paths_tas], ["tas"]; filename_format)

# load different variables from same model (TODO: dimension should be adapble too)
data = mwd.loadPreprocData([paths_tas[end-6], paths_tos[3]], filename_format; dtype="cmip")


# TODO: merge same models (on exisitng dimension) not yet implemented
base = "/albedo/work/projects/p_pool_clim_data/CMIP6/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Amon/pr/gn/v20200212/"
paths = [
    joinpath(base, "pr_Amon_AWI-ESM-1-1-LR_historical_r1i1p1f1_gn_185001-185012.nc"),
    joinpath(base, "pr_Amon_AWI-ESM-1-1-LR_historical_r1i1p1f1_gn_185101-185112.nc")
]
# data = mw.loadPreprocData(paths, :cmip)