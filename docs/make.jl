using Documenter, ModelWeights

DocMeta.setdocmeta!(ModelWeights, :DocTestSetup, :(using ModelWeights); recursive=true)

makedocs(
    sitename = "ModelWeights.jl",
    modules = [ModelWeights],
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "Examples" => "examples.md",
        "weights.md",
        "Manual" => ["Manual/loading-data.md", "Manual/filtering-data.md", "Manual/computing-weights.md"],
        #"Examples" => ["examples/climwip.md", "examples/lgm.md"],
        "Development" => "dev-guide.md",
        "References" => "references.md"
    ],
    format = Documenter.HTML(
        edit_link = "main",
        assets = ["assets/custom.css"],
    ),
    clean = true
)

deploydocs(
    repo = "github.com/awi-esc/ModelWeights.jl.git",
    devbranch = "main"
)
