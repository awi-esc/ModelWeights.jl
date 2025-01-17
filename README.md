# ModelWeights

This repository contains the code for the Julia Package ModelWeights, which computes so called **independence** and **performance** weights for each loaded climate model, namely according to the approach from Brunner et al. (2020). 
The documentation is found here: 

<!-- This approach assigns lower weights to models that make similar predictions (in a predefined reference period). Models whose predictions are further away (compared to the predictions of the other models in the ensemble) receive larger weights.  -->

### Getting started: Initialize project
To activate the Julia project, do ```using Pkg; Pkg.activate(".");``` in the top-level directory or 
alternatively ```pkg> activate .``` (from the package mode, entered in Julia REPL by typing ']'). 
This just activates the project, i.e. it makes it the active project.

Instead of running ```pkg> activate```, you can also specify the project on startup using --project='', e.g. for running a specific script:
```
julia --project=. scripts/plotting/calculate_weights_climwip.jl
```

Then run```pkg> instantiate```. When there is no Manifest.toml file in the project, this will install the latest versions of the dependencies compatible with the project and a Manifest.toml file will be generated (I think this can also be done explicitly by running ```pkg> resolve```). If the project does have a Manifest.toml file, it will install the packages in exactly the same state as given by the manifest file.

To update the packages specified in the Manifest.toml file, run ```pkg> update``` and then run ```instantiate``` again. The Manifest.toml file contains the resolved versions of all dependencies.

Note: ```pkg> resolve```  *updates the dependency graph*n and potentially updates the Manifest.toml file but does not necessarily *install* new compatible versions of dependencies. Contrary to that, ```pkg> update``` fetches and installs the latest versions of packages compatible with constraints in the Project.toml file. 



###  Structure

- scripts/: contains Julia scripts, to show examples and for trying out

- src/: directory contains source code, including all explicitly exported functions for the module

- test/: directory contains tests

- reproduce-climwip-figs
    - recipe_climwip_test_basic_data: contains the output from EsmValTool when running the recipe recipe_climwip_test_basic.
    The directory 'work' contains the processed data that is used by ESMValTool for the diagnostics.
    - plots-replicated-with-julia/: directory that contains the reproduced figures
    - scripts:
        - plotting: data from the work-directory, which is output from ESMValTool and which contains the processed data, is loaded and used to reproduce the Figures from the climwip recipes available in ESMValTool.
        - compute-data-for-diagnostics-from-preproc.jl:  the data that was output from ESMValTool after the preprocessing (so in output folder 'preproc'), is used to reproduce the data that is eventually used (i.e. the data in the work-directory) for creating the plots. 
    - src: contains Julia functions



### Difference paper vs. ESMValTool implementation

|                    | ESMValTool - test_basic                           | ESMValTool - brunner_20esd                                                        |
|--------------------|---------------------------------------------------|-----------------------------------------------------------------------------------|
|performance weights | Climatology (tasCLIM, prCLIM, pslCLIM)            | Anomaly (tasANOM, pslANOM), Standard deviation (tasSTD, pslSTD), Trend (tasTREND) |
|independence weights| Climatology; tas, pr                              | Climatology; tas, psl       |
|observational data  | ERA5                                              |  ERA5                       |
|models              | subset of CMIP6 with constraints                  | all CMIP6 with constraints  |
|mask                | masks out sea and applied to specific region only | nothing, applied everywhere |


In the recipe_climwip_brunner20esd.yml, two models are excluded that were used in the paper (CAMS-CSM1-0 and MPI-ESM1-2-HR (r2)) and only the ERA5 observerational data is used (see [documentation](https://docs.esmvaltool.org/en/latest/recipes/recipe_climwip.html#brunner-et-al-2020-recipe-and-example-independence-weighting)).

Constraints for CMIP6 models: They must provide surface air temperature (tas), and sea level pressure (psl) for the historical, SSP1-2.6 and SSP5-8.5-experiments. 

In the paper,
- they use models post-processed within ETH Zurich CMIP6 next generation archive (provides additional quality checks).
- to represent historical observations in tas and psl, they use two reanalysis products: ERA5 and MERRA-2, each regridded using 2nd order conservative remapping, evaluated in period from 1980-2014. More specifically, they use the mean of both for each grid point.
- they apply the weighting to projections of annual-mean global-mean temperature change from two SSPs (weak, SSP1-2.6 and strong, SSP5-8.5): periods 2041-2060 and 2081-2100 are compared to baseline of the period from 1995-2014
- weights are calculated based on annual mean data from 1980-2014


## Tests
To run the testsets, make sure you activate the Project of your package (e.g. ModelWeights), then go to the pkg> mode (by typing ']' in the REPL) and run 'test'. This starts the tests loaded in 'test/runtests.jl'.

## References
- Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner, Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming from CMIP6 Projections When Weighting Models by Performance and Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020): 995–1012. https://doi.org/10.5194/esd-11-995-2020.

- For online documentation of climwip in ESMValTool, see [here](https://docs.esmvaltool.org/en/latest/recipes/recipe_climwip.html).

## Conceptual thoughts

#### Independence and performance weights
In the ESMVAlTool climwip recipes documentation, they say: 

==Warning: Using only the independence weighting without any performance weighting might not always lead to meaningful results! The independence weighting is based on model output, which means that if a model is very different from all other models as well as the observations it will get a very high independence weight (and also total weight in absence of any performance weighting). This might not reflect the actual independence. It is therefore recommended to use weights based on both independence and performance for most cases.==

