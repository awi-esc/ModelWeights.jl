# ModelWeights

This repository contains the code for the Julia Package ModelWeights. Besides functionalities to load, subset and manipulate CMIP data, it provides functions to compute **independence** and **performance** weights according to the approach from Brunner et al. (2020). 

<!-- This approach assigns lower weights to models that make similar predictions (in a predefined reference period). Models whose predictions are further away (compared to the predictions of the other models in the ensemble) receive larger weights.  -->

The documentation is [here](https://awi-esc.github.io/ModelWeights.jl/).


### Getting started: Initialize project
To activate the Julia project, do ```using Pkg; Pkg.activate(".");``` in the top-level directory or 
alternatively ```pkg> activate .``` (from the package mode, entered in Julia REPL by typing ']'). 
This just activates the project, i.e. it makes it the active project.

Instead of running ```pkg> activate```, you can also specify the project on startup using --project='', e.g. for running a specific script:
```
julia --project=. path_to_your_julia_script.jl
```

Then run```pkg> instantiate```. When there is no Manifest.toml file in the project, this will install the latest versions of the dependencies compatible with the project and a Manifest.toml file will be generated (I think this can also be done explicitly by running ```pkg> resolve```). If the project does have a Manifest.toml file, it will install the packages in exactly the same state as given by the manifest file.

To update the packages specified in the Manifest.toml file, run ```pkg> update``` and then run ```instantiate``` again. The Manifest.toml file contains the resolved versions of all dependencies.

Note: ```pkg> resolve```  *updates the dependency graph*n and potentially updates the Manifest.toml file but does not necessarily *install* new compatible versions of dependencies. Contrary to that, ```pkg> update``` fetches and installs the latest versions of packages compatible with constraints in the Project.toml file. 



###  Structure
- scripts/: directory that contains Julia scripts, to show examples and for trying out

- src/: directory that contains julia source code

- test/: directory that contains tests

- reproduce-climwip-figs
    - recipe_climwip_test_basic_data: contains the output from EsmValTool when running the recipe recipe_climwip_test_basic.
    The directory 'work' contains the processed data that is used by ESMValTool for the diagnostics.
    - plots-replicated-with-julia/: directory that contains the reproduced figures
    - scripts:
        - plotting: data from the work-directory, which is output from ESMValTool and which contains the processed data, is loaded and used to reproduce the Figures from the climwip recipes available in ESMValTool.
        - compute-data-for-diagnostics-from-preproc.jl:  the data that was output from ESMValTool after the preprocessing (so in output folder 'preproc'), is used to reproduce the data that is eventually used (i.e. the data in the work-directory) for creating the plots. 
    - src: contains Julia functions


## Tests
To run the testsets, make sure you activate the Project of your package (e.g. ModelWeights), then go to the pkg> mode (by typing ']' in the REPL) and run 'test'. This starts the tests loaded in 'test/runtests.jl'.

## References
- Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner, Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming from CMIP6 Projections When Weighting Models by Performance and Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020): 995–1012. https://doi.org/10.5194/esd-11-995-2020.

- For online documentation of climwip in ESMValTool, see [here](https://docs.esmvaltool.org/en/latest/recipes/recipe_climwip.html).
