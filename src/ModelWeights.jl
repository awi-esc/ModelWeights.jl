module ModelWeights

# include submodules
include("Data.jl")
include("Timeseries.jl")
include("Weights.jl")
include("Plots.jl")

# make submodules available in the scope of main module (ModelWeights)
import .Data  
import .Timeseries
import .Weights
import .Plots

# export statements: exported functions can be used without qualified name after running 'using ModelWeights'.


end # module ModelWeights
