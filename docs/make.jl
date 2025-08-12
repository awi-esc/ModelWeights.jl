push!(LOAD_PATH,"../src/")

using Documenter, ModelWeights

DocMeta.setdocmeta!(ModelWeights, :DocTestSetup, :(using ModelWeights); recursive=true)

makedocs(
    sitename = "ModelWeights.jl",
    modules = [ModelWeights],
    pages = [
        "Home" => "index.md",
        "getting-started.md",
        "weights.md",
        "Manual" => ["Manual/loading-data.md", "Manual/filtering-data.md", "Manual/computing-weights.md"],
        #"Examples" => ["examples/climwip.md", "examples/lgm.md"],
        "Reference" => ["reference/API.md", "reference/references.md"] 
    ],
    format = Documenter.HTML(
        edit_link = "main"
    ),
    clean = true
)

deploydocs(
    repo = "github.com/awi-esc/ModelWeights.jl.git",
    devbranch = "main"
)
