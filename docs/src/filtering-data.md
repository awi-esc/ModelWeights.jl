# Filtering data

## General options


## Additional options for data preprocessed with ESMValTool

### Data structure
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
│       └── possibly other output from ESMValTool
....
```

In ESMValTool, `BASE_DIR` is the name of the recipe plus the date and time using the format YYYYMMDD_HHMMSS.
`ALIAS` corresponds to the name (chosen by the user) of the diagnostics defined in the recipe.
`VAR_STATISTIC` refers to the name of the entries among the `variables` section in the ESMValTool recipe.
We set these names to the short name of the variable (e.g. 'tas') and the applied statistic (e.g. 'CLIM'), 
concatenated by an underscore (e.g, tas_CLIM). Any name can be chosen here; but when using VAR_STAT,
the data can be subset to certain statistics only (see below).


The structure is basically the same if there is not a seperate subdirectory for each climate 
variable, except that the BASE_DIR refers to the directory that immediately contains the 
preproc-subdirectory:



If there is one subdirectory for each climate variable (this is, for instance, the case when one ESMValTool recipe is used for a single variable), the structure of the data is nearly the same. The only difference is that the upper most directory, BASE_DIR, has subdirectories which then point to the directory called 'preproc'. The name of each of the subdirectories is assumened to contain _VAR_, e.g. (BASE_DIR/recipe_tas_lgm or BASE_DIR/recipe_tos_lgm).
If ESMValTool is used, the name of the subdirectories (`subdir1`, `subdir2`) is the name of the recipe concatenated with the timestamp in the format YYYYMMDD_HHMMSS by an underscore.

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
