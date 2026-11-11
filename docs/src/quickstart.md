# Quickstart: Examples

The paths we use in the examples point to the sample test data, stored in `test/data`. This data was preprocessed with ESMValTool.

First, we need to load `ModelWeights.Data` into our julia environment.

```julia
import ModelWeights.Data as mwd
```

For simplicity, we also do (assuming the examples are run from the top-level directory):

```julia
base = "test/data"
```

## Quickstart: Loading Data

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
julia> dm_models = mwd.subsetModelData(dm, :model)
::DataMap
lgm_tas_CLIM: (72, 36, 15)
lgm_tos_CLIM: (72, 36, 16)
```

Subset data so that every dataset contains the exact same model simulations:

```julia-repl
julia> dm_members = mwd.subsetModelData(dm, :member)
::DataMap
lgm_tas_CLIM: (72, 36, 15)
lgm_tos_CLIM: (72, 36, 15)
```

## Quickstart: Loading + Filtering Data

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
````

To see which models/members are shared across datasets:

```julia-repl
mwd.sharedModels(dm, :model)
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

## Quickstart: Manipulating data

