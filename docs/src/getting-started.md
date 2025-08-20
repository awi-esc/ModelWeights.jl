# Getting started

## Installation

For now, ModelWeights.jl should be used like any local package under development:

1. Clone the repository to a local directory.

2. Add ModelWeights.jl to your workspace in Julia using `Pkg.dev()` rather than `Pkg.add()`:

    ```julia
    using Pkg; Pkg.dev("path/to/ModelWeights.jl")
    ```

3. Now use it like a normal package. If you will make changes to it, then consider also
loading Revise.jl (and do so before using ModelWeights.jl):

    ```julia
    using Revise
    using ModelWeights
    ```


## Basics
