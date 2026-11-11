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





## Documentation

- stable: [https://awi-esc.github.io/ModelWeights.jl/stable/](https://awi-esc.github.io/ModelWeights.jl/stable/)
- dev: [https://awi-esc.github.io/ModelWeights.jl/dev/](https://awi-esc.github.io/ModelWeights.jl/dev/)



## References
- Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner, Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming from CMIP6 Projections When Weighting Models by Performance and Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020): 995–1012. https://doi.org/10.5194/esd-11-995-2020.

- For online documentation of climwip in ESMValTool, see [here](https://docs.esmvaltool.org/en/latest/recipes/recipe_climwip.html).
