push!(LOAD_PATH,"../src/")

using Documenter, ModelWeights

DocMeta.setdocmeta!(ModelWeights, :DocTestSetup, :(using ModelWeights); recursive=true)

makedocs(
    sitename = "ModelWeights.jl",
    modules = [ModelWeights],
    pages = [
        "Home" => "index.md",
        "requirements.md",
        "getting-started.md",
        "weights.md",
        "Examples" => ["examples/climwip.md", "examples/lgm.md"],
        "references.md"
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
