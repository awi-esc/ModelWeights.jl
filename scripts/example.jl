import SimilarityWeights as sw

base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical";
config_path = "/albedo/home/brgrus001/SimilarityWeights/configs/recipe_configs_historical";

# by default all data from config yml files located at config_path is loaded
data_all = sw.loadData(base_path, config_path; dir_per_var=true);
# show all data:
data_all.ids

# Load a subset of the data only
data = sw.loadData(
    base_path,
    config_path;
    isModelData = true,
    dir_per_var = true, 
    common_models_across_vars = false,
    subset = Dict(
        "variables" => ["tas", "pr"], 
        "aliases" => ["historical1"],
        "statistics" => Vector{String}(), # identical to not setting it
        "models" => Vector{String}(), # identical to not setting it
        "data_type" => Vector{String}() # identical to not setting it; e.g. ERA5, CMIP, CMIP5. 
        # loadPreprocData will only load files containing all(!) given values, 
        # if empty and is_model_data=true, only CMIP data will be loaded, if empty and model_data=false,
        # only ERA5 data will be loaded.
    )
);

# can be done directly if argument set to true in loadData above
joint_data = sw.getCommonModelsAcrossVars(data)


# TODO: add examples for computing weights for this data

# amoc_data = sw.loadData(
#     path_to_config_dir, 
#     base_path; 
#     dir_per_var=true, 
#     constraints=sw.DataConstraint(
#         variables = ["msftmz"]
#     )
# )