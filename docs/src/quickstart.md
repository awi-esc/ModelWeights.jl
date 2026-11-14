# Quickstart

The paths we use in the examples point to the sample test data, stored in `test/data`. This data was preprocessed with ESMValTool.

First, we need to load `ModelWeights.Data` into our julia environment.

```julia
import ModelWeights.Data as mwd
```

For simplicity, we also do (assuming the examples are run from the top-level directory):

```julia
base = "test/data"
```

## Loading Data

Load `tos` and `tas` climatologies for the `lgm` experiment:

```julia
paths_lgm_tas = joinpath.(base, ["lgm-cmip5-tas-climatologies", "lgm-cmip6-tas-climatologies"])
paths_lgm_tos = joinpath.(base, ["lgm-cmip5-tos-climatologies", "lgm-cmip6-tos-climatologies"])
ids = ["lgm_tas_CLIM", "lgm_tos_CLIM"] # your chosen ids to refer to the data
```

```julia-repl
julia> dm = mwd.defineDataMap([paths_lgm_tas, paths_lgm_tos], ids; filename_format = :esmvaltool)
::DataMap
lgm_tas_CLIM: (72, 36, 17)
lgm_tos_CLIM: (72, 36, 16)
```

Subset data so that every dataset contains the same models (but not necessarily the same model members):

```julia-repl
julia> dm_models = mwd.subsetModelData(dm; level = :model)
::DataMap
lgm_tas_CLIM: (72, 36, 15)
lgm_tos_CLIM: (72, 36, 16)
```

Subset data so that every dataset contains the exact same model simulations:

```julia-repl
julia> dm_members = mwd.subsetModelData(dm; level = :member)
::DataMap
lgm_tas_CLIM: (72, 36, 15)
lgm_tos_CLIM: (72, 36, 15)
```

## Loading + Filtering Data

The same filtering, as done in the previous examples, can also be done directly when loading the data: 

```julia-repl
julia> mwd.defineDataMap(
    [paths_lgm_tas, paths_lgm_tos], 
    ids; 
    filename_format = :esmvaltool,
    level = :model
)
::DataMap
lgm_tas_CLIM: (72, 36, 15)
lgm_tos_CLIM: (72, 36, 16)
```

for exact simulations, we set `level=:member`:

```julia-repl
julia> mwd.defineDataMap(
    [paths_lgm_tas, paths_lgm_tos], 
    ids; 
    filename_format = :esmvaltool,
    level = :member
)
::DataMap
lgm_tas_CLIM: (72, 36, 15)
lgm_tos_CLIM: (72, 36, 15)
```

To see which models/members are shared across datasets:

```julia-repl
julia> mwd.sharedModels(dm, :model)
13-element Vector{String}:
 "AWI-ESM-1-1-LR"
 "CESM2-FV2"
 "CESM2-WACCM-FV2"
 "CNRM-CM5"
 "FGOALS-g2"
 "GISS-E2-R"
 "INM-CM4-8"
 "IPSL-CM5A-LR"
 "MIROC-ES2L"
 "MIROC-ESM"
 "MPI-ESM-P"
 "MPI-ESM1-2-LR"
 "MRI-CGCM3"
```

```julia-repl
julia> mwd.sharedModels(dm, :member)
15-element Vector{String}:
 "AWI-ESM-1-1-LR#r1i1p1f1_gn"
 "CESM2-FV2#r1i2p2f1_gn"
 "CESM2-WACCM-FV2#r1i2p2f1_gn"
 "CNRM-CM5#r1i1p1"
 "FGOALS-g2#r1i1p1"
 "GISS-E2-R#r1i1p150"
 ⋮
 "MIROC-ES2L#r1i1p1f2_gn"
 "MIROC-ESM#r1i1p1"
 "MPI-ESM-P#r1i1p1"
 "MPI-ESM-P#r1i1p2"
 "MPI-ESM1-2-LR#r1i1p1f1_gn"
 "MRI-CGCM3#r1i1p1"
```

## Manipulating data

We may want to apply a function to every entry in our `DataMap`. For example, to compute (area weighted) global means for all datasets:

```julia 
gms = mwd.apply(dm_members, mwd.globalMeans)
::DataMap
lgm_tas_CLIM: (15,)
lgm_tos_CLIM: (15,)
```

Tpply a function just on a subset of the datasets in our `DataMap`:

```julia-repl
mwd.apply(
    dm_members,
    mwd.globalMeans;
    ids = ["lgm_tas_CLIM"]
)
::DataMap
lgm_tas_CLIM: (15,)
```

To save the result as a new entry, use the data-modifying function `apply!`. Note that, if no new id is given, the old id will be overwritten.

```julia-repl
julia> mwd.apply!(
    dm_members, 
    mwd.globalMeans; 
    ids = ["lgm_tas_CLIM"],
    ids_new = ["lgm_tas_GM"]
)

julia> dm_members
::DataMap
lgm_tas_CLIM: (72, 36, 15)
lgm_tas_GM: (15,)
lgm_tos_CLIM: (72, 36, 15)
```

## Computing Weights

!!! info "Coming soon"
    This section is not finished yet. 

Let's first create some random (non realistic) observational data (with the same grid as in the model data!) and load it into a DataMap:

```julia
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
dm_obs = mwd.defineDataMap([obs_tas, obs_tos], ["tas_CLIM", "tos_CLIM"])
```

We need to specify the diagnostics that the performance and independence weights should be based on as well as the relative impact of each and the hyperparameters $\sigma_D, \sigma_S$ that respectively regulate the strength of the performance- and independence weighting.

Note that the diagnostics for the distances between models and data must all be defined on a lon x lat grid.

```julia
config = mww.ConfigWeights(
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
```

```@raw html
<!-- Let's compute the diagnostics we haven't computed yet that we just specified, the anomalies with respect to the global mean (denoted as "ANOM-GM"). For model data:

```julia
dm_models = mwd.subsetModelData(dm, :model)
ids = collect(keys(dm_models))
ids_new = map(x -> replace(x, "CLIM" => "ANOM-GM"), ids)
mwd.apply(dm_models, mwd.anomaliesGM; ids, ids_new)
```

And for observational data:

```julia
ids = collect(keys(dm_obs))
ids_new = map(x -> replace(x, "CLIM" => "ANOM-GM"), ids)
mwd.apply!(dm_obs, mwd.anomaliesGM; ids, ids_new)
```

Before we can compute weights, we need to make sure that the keys for the observations and models for the diagnostics are identical:

```julia
mwd.renameDict!(
    dm_models,
    ["lgm_tas_ANOM-GM", "lgm_tos_ANOM-GM"],
    ["tas_ANOM-GM", "tos_ANOM-GM"]
)
```

Now, we're ready to compute model weights:

```julia
weights = mww.climwipWeights(dm_models, dm_obs, config)
``` -->
```
