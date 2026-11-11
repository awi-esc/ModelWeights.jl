## Tests
To run the testsets, make sure you the Project of the package is activated, then go to the pkg> mode (by typing ']' in the REPL) and run 'test'. This starts the tests loaded in 'test/runtests.jl'.


## Documentation

To directly see changes you make to the documentation, use LiveServer.jl (add it to your default julia environment which is shared by all other environments). The following commands starts a webserver on localhost:8000 where you see your changes: 

```shell
 julia --project=docs -e 'using ModelWeights, LiveServer; servedocs(port=8000, include_dirs=["src"])'
```

##  Project structure

```@raw html
<!-- - scripts/: directory that contains Julia scripts, to show examples and for trying out -->
```

- src/: directory that contains the julia source code of the package
- test/: directory for tests
- docs/: directory for the documentation (using `Documenter.jl`)
- ext/: extension for the plotting module of ModelWeights

```@raw html
<!-- - reproduce-climwip-figs
    - recipe_climwip_test_basic_data: contains the output from EsmValTool when running the recipe recipe_climwip_test_basic.
    The directory 'work' contains the processed data that is used by ESMValTool for the diagnostics.
    - plots-replicated-with-julia/: directory that contains the reproduced figures
    - scripts:
        - plotting: data from the work-directory, which is output from ESMValTool and which contains the processed data, is loaded and used to reproduce the Figures from the climwip recipes available in ESMValTool.
        - compute-data-for-diagnostics-from-preproc.jl:  the data that was output from ESMValTool after the preprocessing (so in output folder 'preproc'), is used to reproduce the data that is eventually used (i.e. the data in the work-directory) for creating the plots. 
    - src: contains Julia functions -->
```


