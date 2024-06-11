# Climwip-Julia

We aim to reproduce the model-weighting approach from Brunner et al. (2020), starting with how it is implemented in ESMValTool, which is a simplified version of the approach/data presented in their paper (for differences see below).
In scripts/plotting, the data from the work-directory, which is output from ESMValTool and which contains the processed data, is loaded and used to reproduce their Figures in Julia. 

In scripts/compute-data-for-diagnostics-from-preproc.jl, the data that was output from ESMValTool after the preprocessing (so in output folder 'preproc'), is used to reproduce the data that is eventually used (i.e. the data in the work-directory) for creating the plots. 

## Getting started
To activate the climwipJulia project, do ```using Pkg; Pkg.activate(".");``` in the top-level directory. 
This just activates the ClimWipJulia project, i.e. it makes it the active project.

If the project contains a Manifest.toml file, running ```instantiate``` (in package manager command line) after activating the project will install the packages in the same state that is given by the manifest file. Otherwise, it will resolve the latest versions of the dependencies compatible with the project. Then, the dependencies will be installed.  
To update the packages specified in the Manifest.toml file, run ```update``` (in package manager command line) and then run ```instantiate``` again.

There is also Pkg.resolve(). This *updates* the dependencies, i.e. if later versions that are compatible with the other packages and versions are available, these will be loaded, but not necessarily installed, this is done by 'instantiate'. And 'resolve' *updates the dependency graph* and potentially updates the Manifest.toml file. Contrary to that, 'update' fetches and installs the latest versions of packages compatible with constraints in the Project.toml file. 

On albedo: 
- julia -e 'using Pkg; Pkg.activate("."); Pkg.update(); Pkg.instantiate();'
- julia --project=. scripts/plotting/calculate_weights_climwip.jl

Alternatively, you can go to the julia package-command interface (by typing ']' in the Julia REPL) where you should see the name of the julia environment in parentheses in the beginning of the line (like in a conda environment) and run 'update' and/or 'instantiate'. 

#### Note on using on Albedo
On Albedo, I installed the latest version of Julia with juliaup. Make sure to use this one by setting the Julia: executablePath in the settings.json file to the respective path.


### Difference paper vs. ESMValTool implementation

|                    | ESMValTool               | [paper](https://github.com/lukasbrunner/ClimWIP) |
| ------------------ | ------------------------ | ------------------------------------------------ |
|performance weights | Climatology (tasCLIM)    | tasANOM, tasSTD, pslANOM, pslSTD, tasTREND       | 
|independence weights| tasCLIM, pslCLIM         | tasCLIM, pslCLIM                                 |
|observational data  | ERA5                     |   mean(ERA5, MERRA-2)                            |
|models              | subset of CMIP6 with constraints | all CMIP6 with constraints               |

Constraints for CMIP6 models: They must provide surface air temperature (tas), and sea level pressure (psl) for the historical, SSP1-2.6 and SSP5-8.5-experiments. 

In the paper,
- they use models post-processed within ETH Zurich CMIP6 next generation archive (provides additional quality checks).
- to represent historical observations in tas and psl, they use two reanalysis products: ERA5 and MERRA-2, each regridded using 2nd order conservative remapping, evaluated in period from 1980-2014. More specifically, they use the mean of both for each grid point.
- they apply the weighting to projections of annual-mean global-mean temperature change from two SSPs (weak, SSP1-2.6 and strong, SSP5-8.5): periods 2041-2060 and 2081-2100 are compared to baseline of the period from 1995-2014
- weights are calculated based on annual mean data from 1980-2014


###  Structure
- recipe_climwip_test_basic_data:
contains the output from EsmValtool when running the recipe recipe_climwip_test_basic.

- scripts:
    - plotting: directory that contains Julia scripts that generate the figures 
    - compute-data-for-diagnostics-from-preproc.jl: code that reconstructs the weights, distance matrices etc. from the (preprocessed) data output from ESMValTool

- src: directory contains Julia helper functions

- plots-replicated-with-julia: directory that contains the reproduced figures