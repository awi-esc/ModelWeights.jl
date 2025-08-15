# Filtering data

When loading data with the function `defineDataMap`, there is the option to filter the data in several ways, using the the optional keyword argument `constraint` which is a dictionary.

All keys in `constraint` except for two map to a Vector of Strings. The two exceptions are 
`level_shared` and `timeseries` discussed seperately below. 

## Filtering based on filenames

The following keys can be used to constrain the data that shall be loaded into the DataMap:

````julia

constraint = Dict(
    "variables" => ["tas", "tos"],
    "models" => [],
    "variants" => [],
    "experiments" => [],
    "mips" => [],
    "grids" => [],
    "table_ids" => [],
    "filename" => [],
    "timeranges" => Dict(),
    "level_shared" => # one of: "member", "model", :member, :model
)

````

- `filename`:
- `mips`: 
- `models`: if provided, only these models are loaded. It is possible to include (1)
general models (e.g., "ACCESS-CM2"), which would load all available simulations of the respective
model or (2) specific model simulations. The latter are specified as the `MODELNAME#VARIANT`, e.g., `"ACCESS-CM2#r1i1p1f1"`.  



## `level_shared`: Datasets that share the same members/models

The function that is called to constrain the loaded data by the models/members shared across all 
datasets is `subsetModelData`. Note that this loads the data into memory.


## `timeseries`: 


### Optional parameters for filtering data
For both functions, `loadDataFromESMValToolRecipes` and `loadDataFromYAML`, there is a set of optional parameters in order to constrain the loaded data:

- `preview`: If set to false (default), the data will not be loaded and only the metadata with the information that we added that specifies which data will be loaded is returned.

- `dir_per_var`: If set to true (default), only subdirectories of the base\_path that contain `_VARIABLE` in their name will be searched. 


```@raw html
<!-- `subset_shared`: Can be either `nothing` (default), `ModelWeights.MODEL` or `ModelWeights.MEMBER`. If set to `MODEL`, only data from the same models will be loaded. If, for instance, data from lgm and historical experiments shall be loaded, this configuration will ensure that for both, only the same models are loaded, i.e. those for which both experiments were done. While this considers models, not specific simulations, setting subset\_shared to `MEMBER` would only load models that share the exact same simulations (i.e. the same member\_id abbreviation, e.g. `r1i1p1f1`). -->
```

- `subset` is an optional parameter of type `Constraint` or `Nothing`. A Constraint further constrains the data to be loaded and has the following fields:
   
    - `variables`: The short name of the climate variable, e.g. ["tas", "tos"].
    
    - `statistics`: The statistic/diagnostic that was used when preprocessing the data. The names are generally arbitrary, but need to be identical to those you used in the preprocessing; for instance, we refer to the climatological average as "CLIM", s.t. an example value for this argument is ["CLIM"].
    
    - `aliases`: Like for `statistics`, the names are arbitrary but must be identical to what you used when preprocessing the data; e.g. ["historical", "historical1"].
    An alias should refer to a certain `timerange` of a certain experiment, e.g. we call the time period from 1951-1980 'historical1'. To load only this data, it thus does not matter whether you set `timerange=["1951-1980"]` or `alias=["historical1"]`.
    
    - `timeranges`: Timeranges to be loaded; especially for historical data, you may preprocess data for different timeranges,  e.g. ["full", "1980-2014"]. 

    - `projects`:  e.g. ["CMIP5", "CMIP6"]. All filenames of the data to be loaded must contain at least one of the given strings. If not specified and `is_model_data=true`, it is set to ["CMIP"]. If not specified and `is_model_data=false`, it is set to ["ERA5"], the default observational dataset used.

    - `models`: List of models or individual model members, e.g. ["AWI-CM-1-1-MR"]. All filenames must contain at least one of the given strings + "_". The underscore is important since some models have names that are substrings of other models, e.g. "CNRM-CM5" and "CNRM-CM5-C2". 
   
    - `base_subdirs`: If there is one directory for each climate variable (dir\_per\_var=true), you can use this argument to subset the considered subdirectories: if given, the data will be loaded only from subdirectories of the given base\_dir that contain any of the provided values in `base_subdirs` in their name. That is, it's not necessary to specify full directory names here, parts of it are sufficient.


Note that when using the option to load data preprocessed with ESMValTool from a seperate yaml file, 
in which it's possible to specify filtering options for each datasets, these values will be overwritten,
if the same key is provided in the function argument, too. 
That is, the function argument takes precedence over what had been specified in the yaml file if a key is specified in both.
Otherwise if a key is just given in the yaml file, it will still be applied in the filtering.

