import SimilarityWeights as sw

base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical";
path_to_config_dir = "/albedo/home/brgrus001/SimilarityWeights/configs/recipe_configs";

# by default all data is loaded (as specified in config yml files located at 
# path_to_config_dir)
data_all = sw.loadData(path_to_config_dir, base_path; dir_per_var=true);
# show all data:
data_all.ids

# just load a subset of the data from the configs file located at 'path_to_config_dir'
dc = sw.DataConstraint(
    variables=["tas", "pr"], 
    tasks=["historical1"],
    commonModelsAcrossVars=true
);
data = sw.loadData(path_to_config_dir, base_path; dir_per_var=true , constraints=dc);


# TODO: add examples for computing weights for this data

# amoc_data = sw.loadData(
#     path_to_config_dir, 
#     base_path; 
#     dir_per_var=true, 
#     constraints=sw.DataConstraint(
#         variables = ["msftmz"]
#     )
# )