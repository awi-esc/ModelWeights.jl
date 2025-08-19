# Filtering data

When loading data with the function `defineDataMap`, there is the option to filter the data in several ways, using the the optional keyword argument `constraint`.
Note that the filtering is done based on the information retrieved from the filenames. That is, the provided values must correspond to how the data is stored. For CMIP6 data, variables must, for instance, point to the variables' id, e.g. 'tas', 'psl', etc. 

The following keys are considered for the dictionary `constraint`. If not stated otherwise, they map to a vector.

-  `experiments`: load only data from given experiments, e.g. ["lgm", "historical"]
- `filenames`:  load only data from any of the given files
-  `grids`: load only data with given grids, e.g. ["gn"]
- `mips`: load only data from given projects, e.g. to load only CMIP5 data, set to ["CMIP5"]
- `models`: load only data from given models. It is possible to include (1) general models (e.g., "ACCESS-CM2"), which would load all available simulations of the respective model or (2) specific model simulations. The latter are specified as the `MODELNAME#VARIANT`, e.g., `"ACCESS-CM2#r1i1p1f1"`.  
- `table_ids`: load only data from given tables, e.g. ["Amon", "Omon"]
- `variables`: if provided, only data for given variables is loaded, otherwise all data in the specified directories is loaded.
- `variants`: load only specific runs, e.g. ["r1i1p1f1", ...]
- `timeranges`: load only data from files for given timeranges, e.g. ["185001-1850012", "195001-195012"]
- `level_shared::Union{String, Symbol}`: one of: "member", "model", :member, :model. The function that is called to constrain the loaded data by the models/members shared across all datasets is `subsetModelData`. For more details see next paragraph. 

## `level_shared`: Datasets that share the same members/models

This is the most important constraint. It allows to filter the data such that every loaded dataset contains only data from the same models or, more restricted, the same model variants/simulations. If set to `"member"` (or :member), data is only loaded if the respective \emph{model simulation} is present for all specified datasets, i.e. for every variable and experiment. If it is set to `"model"` (or :model), data is only loaded from models that provide data for all specified datasets, yet independent of the exact simulations. That is, every loaded dataset may have of a different number of simulations for each model.


## Further filtering options for data preprocessed with ESMValTool

Besides the possibility to specify the paths to the data directories directly, there are two more ways to load data that was preprocessed with ESMValTool, using the function `defineDataMap`:

- `defineDataMap(path_data, path_recipes, source; dir_per_var::Bool=true)`: The first argument points to the top-level directory that contains the data, the second argument points to a directory that contains one or more ESMValTool recipes and source must be set to `:esmvaltool_recipes`. See [loading data section based on recipes](@ref loading-data-recipes)  for more information.
The additional keyword argument `dir_per_var` specifies whether there was one recipe per variable. If that's the case, there is a single directory for every variable and thus, to load data for a specific variable only subdirectories of the directory that `path_data` points to that contain the respective variable, more precisely `_VARIABLE` (e.g. `_tas`), in their name will be considered.


- `defineDataMap(yaml_content::Dict)` and `defineDataMap(path_config::String)`: Here, the only positional argument is either a dictionary with the configuration or a path to a yaml file with the respective configurations. See [loading data section based on config files](@ref loading-data-config) for more information.

In both cases, it is possible to constrain the directories from where data is loaded by specifying `base_subdirs` in the dictionary input argument `constraint`:

- `"base_subdirs"::Vector{String}`: If specified and if there is a single directory for every variable (default: true), only data is loaded from subdirectories of the given data path (`base_path_data` + `base_dir` when loading from yaml file) that contain any of the given strings in their name.


TODO: following must be updated
- `subset`
   
    - `variables`: The short name of the climate variable, e.g. ["tas", "tos"].
    
    - `statistics`: The statistic/diagnostic that was used when preprocessing the data. The names are generally arbitrary, but need to be identical to those you used in the recipe for preprocessing the data; for instance, we refer to the climatological average as "CLIM", s.t. an example value for this argument is ["CLIM"].
    
    - `aliases`: Like for `statistics`, the names are arbitrary but must be identical to the name you used in the recipe for preprocessing the data; an alias may, for instance, refer to a certain `timerange` of a certain experiment, e.g. `"historical-1950-1980"`.
    
    - `timeranges`: Timeranges to be loaded; especially for historical data, you may preprocess data for different timeranges,  e.g. ["full", "1980-2014"]. 
   

Note that when using the option to load data preprocessed with ESMValTool from a seperate yaml file, in which it's possible to specify filtering options for each datasets, these values will be overwritten, if the same key is provided in the function argument, too. That is, the function argument takes precedence over what had been specified in the yaml file if a key is specified in both. If a key is just given in the yaml file, it will still be applied in the filtering.

