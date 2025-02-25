module ModelWeights

const MODEL_MEMBER_DELIM = "#"
const MODEL_NAME_FIXES = Dict(
    "FGOALS_g2" => "FGOALS-g2",
    "ACCESS1.3" => "ACCESS1-3"
)

include("data-utils.jl")
include("data-functions.jl")
include("plot-utils.jl")
include("weights.jl")
include("plot-weights.jl")
include("plot-data.jl")
include("main-functions.jl")



export loadDataFromYAML, loadDataFromESMValToolConfigs, computeWeights

end # module ModelWeights
