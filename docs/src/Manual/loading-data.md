# Loading data

`ModelWeights` supports loading CMIP model data as well as observational data. It's main purpose, however, is to facilitate working with ensembles of climate models, thus the focus lies on CMIP data.

So far only data on a longitude-latitude grid is supported. The dimensions of the loaded data are (in this order): `:lon`, `:lat`, `:time` (just for timeseries data) and the last dimension is `:model` for observational data or `:member` for CMIP data where `:member` refers to specific model runs, i.e., variants/simulations.

```@raw html
<!-- Also supported is 3D AMOC data, with a fourth dimension called `:lev` representing depth. -->
```

#### Data representation
The basic representation of data is called a `DataMap` which is a simple dictionary mapping from ids 
(Strings) to the data stored as [YAXArrays](https://juliadatacubes.github.io/YAXArrays.jl/v0.6.1/).

#### Metadata
The metadata of each dataset (the `.properties` Dictionary of each YAXArray in a DataMap) contains the metadata loaded from the original model files.

```@raw html
<!-- as well as additional information if it is provided in the optional keyword argument `meta_data` of function `defineDataMap` (a dictionary with keys of type String). -->
```

## `defineDataMap`: optional arguments

The main entry point to load data, is the function `defineDataMap`, which loads data from specified paths. To load data based on configuration files, there are the functions `loadDataFromESMValToolRecipes` and `loadDataFromYAML` explained further below. 

They all share the following **optional keyword arguments**:

- `constraint`: a `Dict{Symbol, <:AbstractArray{String}}` with certain options to filter the data to be loaded; see [here](filtering-data.md#filtering-data) for the available options.
- `constraint_ts`: a `NamedTuple{(:start_year, :end_year), <:Tuple{Integer, Integer}}` specifying a start and end year for time series data. By default, there are no restrictions, i.e. all time points are included.
- `level`: a `Symbol` out of: `:member, :model, :none`. Specifies on which level the data of the different experiments or variables should agree. If set to `:member`, only the subset of exact model runs that are available for every experiment/variable are loaded. If set to `:model`, the exact runs do not matter, but each loaded model is loaded for every experiment/variable. If set to `:none`, for every experiment/variable all available data is loaded, without any restrictions.
- `is_cmip`: a Boolean that, if `true` (default) assumes certain meta data to be present in the data.
- `filename_format`: a Symbol (default: `:esmvaltool`).
- `sorted`: a Boolean that, if `true` (default) the returned data is sorted alphabetically with respect to the model names.
- `model_times`: a Boolean that, if set to `true` (default: `false`) assumes that the time dimension refers to model times, i.e., not to actual dates. <!--This has an influence on assumptions for loading time series data-->
- `filename_format`: a Symbol that specifies the structure of the filenames of the data (default: `:esmvaltool`). This argument is important since the data is loaded based on the information retrieved from the filenames of the data. The predefined filename formats are explained in the next paragraph.

## Filename formats
  - `filename_format = :esmvaltool` (default)
    Only esmvaltool format is assumed, for every file it is automatically derived which exact format it corresponds to.

  - `filename_format = :esmvaltool_cmip5`
    ```julia
    "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE"
    # Example: CMIP5_ACCESS1-0_Amon_piControl_r1i1p1_tas.nc
    ```


  - `filename_format = :esmvaltool_cmip6`
    ```julia
    "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE_GRID"
    # Example: CMIP6_ACCESS-CM2_Amon_piControl_r1i1p1f1_tas_gn.nc
    ```

  - `filename_format = :esmvaltool_obs`
    ```julia
    "GRID_MODEL_TYPE_VERSION_TABLEID_VARIABLE"
    # Example: 
    ```
    When the data was filtered for a specific time period within ESMValTool, the filename of the stored data ends with the timerange. Such data is loaded with the following filename formats:
  
  - `filename_format = :esmvaltool_cmip5_tr`
    ```julia
    "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE_TIMERANGE"
    # Example: CMIP5_ACCESS1-0_Amon_historical-rcp85_r1i1p1_tas_1960-2100.nc
    ```

  - `filename_format = :esmvaltool_cmip6_tr`
    ```julia
    "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE_GRID_TIMERANGE"
    # Example: CMIP6_ACCESS-CM2_Amon_historical_r1i1p1f1_tas_gn_1850-1900.nc
    ```

  - `filename_format = :esmvaltool_obs_tr`
    ```julia
    "GRID_MODEL_TYPE_VERSION_TABLEID_VARIABLE_TIMERANGE"
    # Example:
    ```


  - `filename_format = :cmip` This is the format of CMIP data as downloaded from ESGF. For every file it is automatically derived which exact cmip format it corresponds to (CMIP5/CMIP6) and whether it is a timeseries or a climatology.

  - `filename_format = :cmip6`
    ```julia
    "VARIABLE_TABLEID_MODEL_EXPERIMENT_VARIANT_GRID"
    # Example: orog_fx_AWI-ESM-1-REcoM_historical_r1i1p1f1_gn.nc
    ```
  - `filename_format = :cmip5`
    ```julia
      "VARIABLE_TABLEID_MODEL_EXPERIMENT_VARIANT"
        # Example:
    ```
  
  - `filename_format = :cmip6_tr`
    ```julia
    "VARIABLE_TABLEID_MODEL_EXPERIMENT_VARIANT_GRID_TIMERANGE"
    # Example: pr_Amon_AWI-ESM-1-1-LR_lgm_r1i1p1f1_gn_396101-396512.nc
    ```
  - `filename_format = :cmip5_tr`
    ```julia
      "VARIABLE_TABLEID_MODEL_EXPERIMENT_VARIANT_TIMERANGE"
        # Example: pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc
    ```

```@raw html
<!-- Further, it is possible to specify your own filename format. TODO  -->
```



## `defineDataMap`: Loading data from directories

When loading data from specified paths, it is generally assumed that the paths point to directories that contain data for a single climate variable and experiment: the dimensions of the loaded data are :lon, :lat, :time (for timeseries) and :model (or :member).

There are four methods of the function `defineDataMap`. The first of two positional arguments is the path (or paths) to where the data is stored and the second argument is the id (or ids) that will be used for the respective dataset(s).

1. For a single dataset, i.e. one experiment and one variable only. For example, to load `tos` model data for the `lgm` (in this example the output is a dataset with 6 models on a 72x36 grid):
   ```julia-repl
   julia> defineDataMap("path/to/lgm/tos", "lgm-tos")
   ::DataMap
   lgm-tos: (72, 36, 6)
   ```
   The first argument should point to the directory that only contains `tos` data for the `lgm`.

2. As in 1., but the data for the requested dataset should be retrieved from multiple locations, e.g. when CMIP5 and CMIP6 are stored separately:
   ```julia-repl
   julia> paths = ["path/to/cmip5/lgm/tos", "path/to/cmip6/lgm/tos"]
   julia> defineDataMap(paths, "lgm-tos")
   ::DataMap
   lgm-tos: (72, 36, 6)
   ```

3. To load data for different datasets, e.g. `tos` and `tas` for the `lgm` experiment:
   ```julia-repl
   julia> paths = ["path/to/lgm/tas", "path/to/lgm/tos"]
   julia> ids = ["lgm_tas", "lgm_tos"]
   julia> defineDataMap(paths, ids)
   ::DataMap
   lgm_tas: (72, 36, 14)
   lgm_tos: (72, 36, 14)
   ```

4. Same as in 3., but the data of one or more datasets is retrieved from multiple locations, e.g., when you store CMIP5 and CMIP6 data separately:
   ```julia-repl
   julia> paths = [
     ["path/to/cmip5/lgm/tas", "path/to/cmip6/lgm/tas"], 
     ["path/to/cmip5/lgm/tos", "path/to/cmip6/lgm/tos"]
   ]
   julia> ids = ["lgm_tas", "lgm_tos"]
   julia> defineDataMap(paths, ids)
   ::DataMap
   lgm_tas: (72, 36, 14)
   lgm_tos: (72, 36, 14)
   ```

## `loadDataFromESMValToolRecipes`: [Loading data from ESMValTool recipes](@id loading-esmvt)

If you use ESMValTool to preprocess data, you can use the recipes for preprocessing the data as config files to load the data (or a subset of it) into a `DataMap`. 

```julia
using ModelWeights
loadDataFromESMValToolRecipes(path_data, path_recipes)
```

The first argument, `path_data`, points to the top-level directory where the data is stored, corresponding to ESMValTool's `output_dir`. The second argument, `path_recipes`, points to a directory that contains one or more ESMValTool recipes. Note that these should correspond to the explicitly spelled out versions (`X_filled.yml`, where X is the name of the original recipe) that are automatically generated by ESMValTool, stored in the `run` directory of `output_dir`.

The following additional *optional* argument can further be used when loading data based on ESMValTool recipes:

- `dir_per_var`: if true (default), every climate variable to be loaded is stored in a subdirectory of `path_data`. This is the case, when each recipe loads data for one variable only.

When the data was preprocessed with ESMValTool, the structure of the directories where the data is stored looks as follows.
Uppercase names refer to variables, lowercase names mean that the names can be arbitrary unless stated differently in the comments to the right. 

```bash
├── BASE_DIR # arbitrary (ESMValTool recipe name + _ + YYYYMMDD_HHMMSS)
│       └── preproc # this directory must be called 'preproc'
│       │      └── ALIAS # arbitrary (name of the diagnostic used in the ESMValTool recipe, chosen by the user)
│       │      │       └── VAR_STATISTIC # name of entry among `variables` in the ESMValTool recipe e.g. tos_CLIM (chosen by the user)
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      │            └── ...
│       │      └── ALIAS # arbitrary (name of the diagnostic used in the ESMValTool recipe, chosen by the user)
│       │      │       └── VAR_STATISTIC # name of entry among `variables` in the ESMValTool recipe e.g. tos_CLIM (chosen by the user)
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      │            └── ...
│       └── other output from ESMValTool
....
```

`BASE_DIR` is the name of the recipe plus the date and time using the format YYYYMMDD_HHMMSS.
`ALIAS` corresponds to the name (chosen by the user) of the diagnostics defined in the recipe.
`VAR_STATISTIC` refers to the name of the entries among the `variables` section in the ESMValTool recipe.
We set these names to the short name of the variable (e.g. 'tas') and the applied statistic (e.g. 'CLIM'), 
concatenated by an underscore (e.g, tas_CLIM). Any name can be chosen here; but when using VAR_STAT,
the data can be subset to certain statistics only (see TODO).


If there is one subdirectory for each climate variable (this is, for instance, the case when one ESMValTool recipe is used for a single variable), the structure of the data is nearly the same. The only difference is that the upper most directory, BASE_DIR, has subdirectories which then point to the directory called 'preproc'. The name of the subdirectories (`subdir1`, `subdir2`) is the name of the recipe concatenated with the timestamp in the format YYYYMMDD_HHMMSS by an underscore.
Note that here, we assume one subdirectory per climate variable, and thus the name of each of the subdirectories is assumened to contain _VAR_, e.g. (BASE_DIR/recipe_tas_lgm or BASE_DIR/recipe_tos_lgm).

```bash
├── BASE_DIR
│   └── subdir1  # name must contain "_VAR_" where VAR is the short_name of the variable, e.g. _tos_
│   └── subdir2  # name must contain "_VAR_" where VAR is the short_name of the variable, e.g. _tos_
│       └── preproc # this directory must be called 'preproc'
│       │      └── ALIAS # the name of the diagnostic in the ESMValTool recipe
│       │      │       └── VAR_STATISTIC  # can be any name, but this schema useful for subsetting easily wrt statistics
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      │            └── ...
│       │      └── ALIAS # the name of the diagnostic in the ESMValTool recipe
│       │      │       └── VAR_STATISTIC  # can be any name, but this schema useful for subsetting easily wrt statistics
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      │            └── ...
│       └── possibly other output from ESMValTool
....
```

### `loadDataFromYAML`: Using .yaml files as config files, including constraints
!!! info "Coming soon"
    This section is under development.


```@raw html
<!-- If the data was preprocessed with ESMValTool, it is also possible to create a config yaml file with the required values. Here is an example that loads `tas` and `tos` data for the experiments `lgm`, `historical` and `piControl`: 

````yaml
base_path_data: "/path/to/dir/with/preprocessed/data"

datasets: [
  {
    base_dir: "LGM", 
    exp: "lgm", 
    variables: ["tas", "tos"],
    aliases: ["lgm"],
    statistics: ["CLIM"], 
    subset: {
      base_subdirs: ["20241114"]
    }
  },
  {
    base_dir: "piControl", 
    exp: "piControl", 
    variables: ["tas", "tos"],
    aliases: ["piControl"],  
    statistics: ["CLIM"],
    subset: {
      base_subdirs: ["lgm_models", "20250322"]
    }
  },
    {
    base_dir: "historical", 
    exp: "historical", 
    variables: ["tas", "tos"],
    aliases: ["historical"], 
    statistics: ["CLIM"]
  }
]
````

Keys expected in the config-yaml file:
- `path_data`: base path to where the preprocessed data is stored. This path will be concatenated with the directories specified for each dataset. It can 
be omitted if inside datasets (see below), `base_dir` refers to the entire absolute path.
- `datasets`: list of each dataset to be loaded. Each dataset is specified by a dictionary with the following keys:
    - `base_dir::String`: relative path to the data with respect to the path given in `path_data`. If `path_data` is not specified in the config file, `path_data` for this dataset is assumed to refer to an absolute path instead.
    - `exp::String`: e.g., "historical"
    - `variables::Vector{String}`: e.g., ["tos", "tas"]
    - `statistics::Vector{String}`: may be empty for fixed variables for which no statistics are computed. In that case, `statistics` is set to ["none"] and a warning is thrown unless for known fixed variables (for now just [`orog`]).

   The following keys are optional:
   - `subdirs`: vector of strings; only subdirectories that contain any of the given strings in their name will be considered.
   - `timeranges`:
   - `aliases`: -->
```
