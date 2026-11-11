# ModelWeights

This repository contains the code for the Julia Package ModelWeights. Besides functionalities to load, subset and manipulate CMIP data, it provides functions to compute **independence** and **performance** weights according to the approach from Brunner et al. (2020).

<!-- This approach assigns lower weights to models that make similar predictions (in a predefined reference period). Models whose predictions are further away (compared to the predictions of the other models in the ensemble) receive larger weights.  -->


## Installation

For now, ModelWeights.jl should be used like any local package under development:

1. Clone the repository to a local directory.

2. Add ModelWeights.jl to your workspace in Julia using `Pkg.dev()` rather than `Pkg.add()`:

```julia
using Pkg
Pkg.dev("path/to/ModelWeights.jl")
```

## Quickstart

The paths we use in the examples point to the sample test data, stored in `test/data`. This data was preprocessed with ESMValTool.

First, we need to load `ModelWeights.Data` into our julia environment.

```julia
import ModelWeights.Data as mwd
```

For simplicity, we also do (assuming the examples are run from the top-level directory):

```julia
base = "test/data"
```

### Loading Data

Load `tos` and `tas` climatologies for the `lgm` experiment:

```julia
paths_lgm_tas = joinpath.(
    base, 
    ["lgm-cmip5-tas-climatologies", "lgm-cmip6-tas-climatologies"]
)
paths_lgm_tos = joinpath.(
    base, 
    ["lgm-cmip5-tos-climatologies", "lgm-cmip6-tos-climatologies"]
)
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

### Loading + Filtering Data

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

### Manipulating data

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


## Documentation

- stable: [https://awi-esc.github.io/ModelWeights.jl/stable/](https://awi-esc.github.io/ModelWeights.jl/stable/)
- dev: [https://awi-esc.github.io/ModelWeights.jl/dev/](https://awi-esc.github.io/ModelWeights.jl/dev/)



## References
- Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner, Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming from CMIP6 Projections When Weighting Models by Performance and Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020): 995–1012. https://doi.org/10.5194/esd-11-995-2020.

- For online documentation of climwip in ESMValTool, see [here](https://docs.esmvaltool.org/en/latest/recipes/recipe_climwip.html).
