# Filtering data

When loading data, the data can be filtered in several ways, explained below. Note that the filtering is done based on the information retrieved from the filenames. By the default, TODO 

<!-- %That is, the provided values must correspond to how the data is stored. For CMIP6 data, variables must, for instance, point to the variables' id, e.g. 'tas', 'psl', etc.  -->

## Possible filtering options
### [Using kwarg `constraint`](@id constraint-options)

- `filenames`:  load only data from any of the given files
- `variables`: if provided, only data for given variables is loaded, otherwise all data in the specified directories is loaded. 
- `table_ids`: load only data from given tables, e.g. `["Amon", "Omon"]`
- `models`: load only data from given models, e.g. `["ACCESS-CM2", "MPI-ESM-P"]` 
-  `experiments`: load only data from given experiments, e.g. ``["lgm", "historical"]``
- `variants`: load only specific runs, e.g. `["r1i1p1f1", ...]`
-  `grids`: load only data with given grids, e.g. `["gn"]`
- `mips`: load only data from given projects, e.g. to load only CMIP5 data, set to ["CMIP5"]
- `members`: load only specific model runs; requires full model name including run (`MODELNAME#VARIANT`), e.g. `["MPI-ESM-P#r1i1p1"]`

### [Using kwarg `constraint_ts` for timeseries constraints](@id constraint-ts-options)

For timeseries data, the start and end year can be specified using the keyword argument `constraint_ts` which is a `NamedTuple{:start_year, :end_year}`. By default, all data is returned.


### Using kwarg `level`: Datasets that share the same members/models

This is probably the most important constraint. It allows to filter the data such that only models are loaded from which data is available for all requested experiments/variables (`level=:model`), independent of the exact simulations available. We may, for instance, want air temperature data of all models for which this variable is available for the `historical` experiment. If we set `level=:model`, the exact simulations do not matter: one model may, for instance, have 2 simulations for the historical experiment, while another one only has a single run for the historical experiment: 

| Model | Experiment | variable | Variant (simulation id)|
|-------|------------|----------|------------------------|
| $m_1$ | historical | tas      | r1i1p1f1               |
| $m_2$ | historical | tas      | r1i1p1f1, r1i1p1f2     |

If the filtering should be more restricted, `level` can be set to `:member` which requires the exact same simulations. So for the simple example from the table above, the output would be restricted to just the `r1i1p1f1` simulation for both models. 

By default, no filtering takes place (`level=:none`).

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

