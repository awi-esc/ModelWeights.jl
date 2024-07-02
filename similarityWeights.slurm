#!/bin/bash

# load software 
module load julia
module load esmvaltool
module load conda

#conda activate /albedo/soft/sw/conda-sw/esmvaltool

#srun esmvaltool run recipe_climwip_test_basic.yml
#srun esmvaltool run recipe_climwip_brunner20esd.yml

julia -e 'using Pkg; Pkg.activate("."); Pkg.update(); Pkg.instantiate();'

srun julia --project=. scripts/plotting/calculate_weights_climwip.jl
srun julia --project=. scripts/plotting/weighted_temperature_graph.jl
srun julia --project=. scripts/plotting/weighted_temperature_map.jl

srun julia --project=. scripts/compute-data-for-diagnostics-from-preproc.jl