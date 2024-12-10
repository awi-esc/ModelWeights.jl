# Requirements

## Preprocessed data
We assume that you have your preprocessed data ready, so that the data from
the different models can be combined, meaning that all models must for instance 
have the same grid.
We do the preprocessing of the data with ESMValTool. A set of recipes to 
download and preprocess some data can be found in our repository [ESMDataPrep](https://github.com/awi-esc/ESMDataPrep).

The structure of the directories where the preprocessed data is stored is
expected by our tool to follow a certain structure. If there is one subdirectory
for each climate variable (this is, for instance, the case when one ESMValTool
recipe is used for a single variable), the structure should be as outline below. 
Uppercase names refer to variables, lowercase names mean that the names can be
arbitrary, unless stated differently in the comments to the right. 

```bash
├── BASE_DIR
│   └── subdir1  # name must contain "_" + the short_name of the variable, e.g. _tos
│   └── subdir2  # name must contain "_" + the short_name of the variable, e.g. _tos
│       └── preproc # this directory must be called 'preproc'
│       │      └── ALIAS # the name of the diagnostic in ESMValTool, e.g. historical
│       │      │       └── VAR_STATISTIC    # e.g. tos_CLIM
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      │            └── ...
│       │      └── ALIAS # the name of the diagnostic in ESMValTool, e.g. historical1
│       │      │       └── VAR_STATISTIC    # e.g. tos_STD
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      │            └── ...
│       └── possibly other output from ESMValTool
....
```

The structure is basically the same if there is not a seperate subdirectory for
each climate variable, except that the BASE_DIR refers to the directory that
immediately contains the preproc-subdirectory: 


```bash
├── BASE_DIR
│       └── preproc # this directory must be called 'preproc'
│       │      └── ALIAS # the name of the diagnostic in ESMValTool, e.g. historical
│       │      │       └── VAR_STATISTIC # e.g. tos_CLIM
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      │            └── ...
│       │      └── ALIAS # the name of the diagnostic in ESMValTool, e.g. historical1
│       │      │       └── VAR_STATISTIC # e.g. tos_STD
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      │            └── ...
│       └── possibly other output from ESMValTool
....
```

## Config files

Further, to load the data, we need one or more yaml configuration files.
From these, we retrieve which combination of variables, statistics, aliases,
experiments and timeranges will be considered.
   
Since we use ESMValTool to preprocess the data, we can simply use our 
ESMValTool recipes (the fully spelled out version output by ESMValTool, which
is stored in the run-folder and has the name of recipe + "_filled.yml") as 
config files here. Note that not everything in the recipes is actually needed 
for loading the data. For a minimal example with the unnecessary sections 
removed, see [this example](https://github.com/awi-esc/SimilarityWeights/blob/main/configs/examples/esmvaltool-recipes/mwe_esmvaltool_config.yml). 

You may have one centralized directory, where you store the all yaml files
(i.e., if you used ESMValTool, all recipes used for preprocessing the data).