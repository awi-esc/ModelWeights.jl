import SimilarityWeights as sw

base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical";
path_to_config_dir = "/albedo/home/brgrus001/SimilarityWeights/configs/recipe_configs";

# specify which data to load (default: all data will be loaded)
dc = sw.DataConstraint(
    variables=["tas"], 
    tasks=["historical1"],
    commonModelsAcrossVars=true
);

data = sw.loadData(path_to_config_dir, base_path; dir_per_var=true , constraints=dc);


# load all data as specified in config yml files
# Here just two variables are used (see folder configs/)
data_all = sw.loadData(path_to_config_dir, base_path; dir_per_var=true);
# show all data:
data_all.ids


# TODO: add examples for computing weights for this data