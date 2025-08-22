# Getting started

## Installation

For now, ModelWeights.jl should be used like any local package under development:

1. Clone the repository to a local directory.

2. Add ModelWeights.jl to your workspace in Julia using `Pkg.dev()` rather than `Pkg.add()`:

    ````julia
    using Pkg; Pkg.dev("path/to/ModelWeights.jl")
    ````

3. Now use it like a normal package. If you will make changes to it, then consider also
loading Revise.jl (and do so before using ModelWeights.jl):

    ````julia
    using Revise
    using ModelWeights
    ````


## Basics

### Quickstart Loading/Filtering Data
Load different datasets into a DataMap

````julia
import ModelWeights as mw

paths_ds1 = ["test/data/lgm-cmip6-tas-climatologies", "test/data/lgm-cmip5-tas-climatologies"]
paths_ds2 = ["test/data/lgm-cmip6-tos-climatologies"]
ids = ["lgm_tas_CLIM", "lgm_tos_CLIM"]
dm = mw.defineDataMap([paths_ds1, paths_ds2], ids; filename_format = :esmvaltool)
````

Subset data so that every dataset contains the same models (but not necessarily the same model members):

````julia
dm_models = mw.subsetModelData(dm, :model)
````

Subset data so that every dataset contains the same model simulations:

````julia
dm_members = mw.subsetModelData(dm, :member)
````

The same filtering can be done directly when loading the data: 

````julia
constraint = Dict(
    "level_shared" => :member
)
dm_members2 = mw.defineDataMap([paths_ds1, paths_ds2], ids; constraint, filename_format = :esmvaltool)
````

To get the shared models/members:

````julia
models_shared = mw.sharedModels(dm, "model")
members_shared = mw.sharedModels(dm, "member")
````

### Quickstart Manipulating Data

To conveniently apply any function that takes as first argument a YAXArray, to every (or a subset) of datasets in a DataMap, you can use the functions `apply` or `apply!`. While `apply` returns a new DataMap, `apply!` modifies the input DataMap.
Both provide keyword arguments `ids` and `ids_new`. If `ids` is not provided, the function is applied to every dataset in the input DataMap and if `ids_new` is not provided, the new names will be identical to the old names. That is, if they are not specified and `apply!` is used the respective entries are overwritten.

For example, to compute (area weighted) global means for all datasets:

```@example 
mw.Data.apply(dm_members, mw.Data.globalMeans)
```

for a subset, modifying original DataMap:

```@example 
mw.Data.apply!(dm_members, mw.Data.globalMeans; ids = ["lgm_tas_CLIM"], ids_new = ["lgm_tas_GM"])
```

The first argument to `apply` (or `apply!`) is always the DataMap, the second the function that takes as first argument a YAXArray, followed by the positional and then possibly keyword arguments that the function expects.


### Quickstart Weights

Let's first create some random (non realistic) observational data (with the same grid as in the model data!) and load it into a DataMap:

````julia
using DimensionalData
using YAXArrays

longitudes = collect(lookup(dm_models["lgm_tas_CLIM"], :lon))
latitudes = collect(lookup(dm_models["lgm_tas_CLIM"], :lat))

obs_tas = YAXArray((
    Dim{:lon}(longitudes), Dim{:lat}(latitudes), Dim{:model}(["obs"])),
    rand(length(longitudes), length(latitudes), 1)
)
obs_tos = YAXArray((
    Dim{:lon}(longitudes), Dim{:lat}(latitudes), Dim{:model}(["obs"])),
    rand(length(longitudes), length(latitudes), 1)
)
dm_obs = mw.defineDataMap([obs_tas, obs_tos], ["tas_CLIM", "tos_CLIM"])
````

We need to specify the diagnostics that the performance and independence weights should be based on as well as the relative impact of each and the hyperparameters $\sigma_D, \sigma_S$ that respectively regulate the strength of the performance- and independence weighting.

Note that the diagnostics for the distances between models and data must all be defined on a lon x lat grid.

````julia
config = mw.Weights.ConfigWeights(
    performance = Dict(
        "tas_ANOM-GM" => 1,
        "tos_ANOM-GM" => 1
    ),
    independence = Dict(
        "lgm_tas_CLIM" => 0.5,
        "lgm_tos_CLIM" => 0.5
    ),
    sigma_independence = 0.5, # may be ommited (0.5 is default value)
    sigma_performance = 0.5 # may be ommited (0.5 is default value)
);
````

Let's compute the diagnostics we haven't computed yet that we just specified, the anomalies with respect to the global mean (denoted as "ANOM-GM"). For model data:

```@example
ids = collect(keys(dm_models))
ids_new = map(x -> replace(x, "CLIM" => "ANOM-GM"), ids)
mw.Data.apply!(dm_models, mw.Data.anomaliesGM; ids, ids_new)
```

And for observational data:

```@example
ids = collect(keys(dm_obs))
ids_new = map(x -> replace(x, "CLIM" => "ANOM-GM"), ids)
mw.Data.apply!(dm_obs, mw.Data.anomaliesGM; ids, ids_new)
```

Before we can compute weights, we need to make sure that the keys for the observations and models for the diagnostics are identical:

```@example
mw.Data.renameDict!(
    dm_models,
    ["lgm_tas_ANOM-GM", "lgm_tos_ANOM-GM"],
    ["tas_ANOM-GM", "tos_ANOM-GM"]
)
```

Now, we're ready to compute model weights:

```@example
weights = mw.Weights.climwipWeights(dm_models, dm_obs, config)
```






