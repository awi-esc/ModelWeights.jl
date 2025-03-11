# Weights for models of lgm experiment


First, we want to load all CMIP5 and CMIP6 model data for the lgm-experiment.

```julia
# points to directory where preprocessed data is stored
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
# point to directory where config files are stored
path_recipes = "/albedo/home/brgrus001/SimilarityWeights/configs/lgm-cmip5-cmip6";

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
model_members_lgm = Array(dims(first(values(lgm_data)), :member))
```

Then, we load data from the same models for the historical experiment:

```julia
base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/";
config_path = "/albedo/home/brgrus001/SimilarityWeights/configs/historical";
historical_data = mw.loadDataFromESMValToolConfigs(
    base_path, config_path;
    only_shared_models = true,
    subset = mw.Constraint(
        statistics = ["CLIM"],
        variables = ["tas", "tos"],
        aliases = ["historical"],
        timeranges = ["full"],
        models = model_members_lgm,
        subdirs = ["20241121", "20241118"]
    )
);
````

This is the resulting data:
```julia
```

We could also have used a single yaml configuration file independent of the ESMValTool recipes and load the data all together one:


```julia
# Load model data for experiment lgm and historical in one run from new config file
path_config = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/configs-SimilarityWeights/example-lgm-historical.yml";

model_data = mw.loadDataFromYAML(
    path_config;
    dir_per_var = true, # default: true
    is_model_data = true, # default: true
    only_shared_models = true, # default: false
    preview = false
);
```


# TODO: show that we get the same outfit different approaches!