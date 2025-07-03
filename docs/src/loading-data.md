# Loading data

## Loading data from directories

## Loading data preprocessed with ESMValTool 

### Based on recipes

### Based on config files

If the data is loaded from ESMValToolRecipes, you may have one centralized directory, where you store the all recipe files, namely those output by ESMValTool (inside the run folder, then 'RECIPE_NAME_filled.yml'). Note that not everything in the recipes is actually needed for loading the data. For a minimal example with the unnecessary sections removed, see [this example](https://github.com/awi-esc/SimilarityWeights/blob/main/configs/examples/esmvaltool-recipes/mwe_esmvaltool_config.yml)


Further, to load the data, we need one or more yaml configuration files. From these, we retrieve which combination of variables, statistics, aliases, experiments and timeranges will be considered.
