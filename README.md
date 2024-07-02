# SimilarityWeights

This repository contains the code for the Julia Module SimilarityWeights, which computes so called *independence weights* for each loaded model, namely according to the approach from Brunner et al. (2020). 
This approach assigns lower weights to models that make similar predictions (in a predefined reference period). Models whose predictions are further away (compared to the predictions of the other models in the ensemble) receive larger weights. 


### Getting started
To activate the Julia project, do ```using Pkg; Pkg.activate(".");``` in the top-level directory. 
This just activates the project, i.e. it makes it the active project.

If the project contains a Manifest.toml file, running ```instantiate``` (in package manager command line) after activating the project will install the packages in the same state that is given by the manifest file. Otherwise, it will resolve the latest versions of the dependencies compatible with the project. Then, the dependencies will be installed.  
To update the packages specified in the Manifest.toml file, run ```update``` (in package manager command line) and then run ```instantiate``` again.

There is also Pkg.resolve(). This *updates* the dependencies, i.e. if later versions that are compatible with the other packages and versions are available, these will be loaded, but not necessarily installed, this is done by 'instantiate'. And 'resolve' *updates the dependency graph* and potentially updates the Manifest.toml file. Contrary to that, 'update' fetches and installs the latest versions of packages compatible with constraints in the Project.toml file. 

#### Note on using on Albedo

On Albedo, I installed the latest version of Julia with juliaup. Make sure to use this one by setting the Julia: executablePath in the settings.json file to the respective path.

Activate and instantiate project by running the following commands from the julia REPL:
```
using Pkg; 
Pkg.activate("."); 
Pkg.update(); 
Pkg.instantiate();'
```
Alternatively, you can go to the julia package-command interface (by typing ']' in the Julia REPL) where you should see the name of the julia environment in parentheses in the beginning of the line (like in a conda environment) and run 'update' and/or 'instantiate'. 


Instead of using 'activate', you can specify the project on startup using --project='', e.g. for running a specific script:
```
julia --project=. scripts/plotting/calculate_weights_climwip.jl
```

###  Structure

- recipe_climwip_test_basic_data/: directory contains the output from ESMValTool for the recipe 'recipe_climwip_test_basic.yml'. The directory 'work' contains the processed data that is used by ESMValTool for the diagnostics.

- scripts/: directory contains scripts for reproducing the figures from the original climwip recipes (inside scripts/plotting). Scripts/compute-data-for-diagnostics-from-preproc.jl loads the preprocessed data (output from ESMValTool, inside preproc-dir) and recomputes the data used for the diagnostics

- src/: 
In scripts/plotting, the data from the work-directory, which is output from ESMValTool and which contains the processed data, is loaded and used to reproduce the Figures from the climwip recipes available in ESMValTool in Julia. 
In scripts/compute-data-for-diagnostics-from-preproc.jl, the data that was output from ESMValTool after the preprocessing (so in output folder 'preproc'), is used to reproduce the data that is eventually used (i.e. the data in the work-directory) for creating the plots. 


- recipe_climwip_test_basic_data:
contains the output from EsmValTool when running the recipe recipe_climwip_test_basic.

- scripts:
    - plotting: directory that contains Julia scripts that generate the figures 
    - compute-data-for-diagnostics-from-preproc.jl: code that reconstructs the weights, distance matrices etc. from the (preprocessed) data output from ESMValTool

- src: directory contains helper functions and in compute-weights.jl the exported functions for the module

- plots-replicated-with-julia/: directory that contains the reproduced figures


- tests: TODO: add tests for exported functions

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


## References
- Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner, Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming from CMIP6 Projections When Weighting Models by Performance and Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020): 995–1012. https://doi.org/10.5194/esd-11-995-2020.

- For online documentation of climwip in ESMValTool, see [here](https://docs.esmvaltool.org/en/latest/recipes/recipe_climwip.html).

## Conceptual thoughts

#### Independence and performance weights
In the ESMVAlTool climwip recipes documentation, they say: 

==Warning: Using only the independence weighting without any performance weighting might not always lead to meaningful results! The independence weighting is based on model output, which means that if a model is very different from all other models as well as the observations it will get a very high independence weight (and also total weight in absence of any performance weighting). This might not reflect the actual independence. It is therefore recommended to use weights based on both independence and performance for most cases.==

