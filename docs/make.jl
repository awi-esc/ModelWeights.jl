using Documenter, DocumenterCitations, ModelWeights

DocMeta.setdocmeta!(ModelWeights, :DocTestSetup, :(using ModelWeights); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs", "refs.bib"); style=:numeric)

makedocs(
    sitename = "ModelWeights.jl",
    modules = [ModelWeights],
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Quickstart" => "quickstart.md",
        "weights.md",
        "Manual" => ["Manual/loading-data.md", "Manual/filtering-data.md", "Manual/computing-weights.md"],
        #"Examples" => ["examples/climwip.md", "examples/lgm.md"],
        "Development" => "dev-guide.md",
        "References" => ["refs/public-api.md", "refs/references.md"]
    ],
    format = Documenter.HTML(
        edit_link = "main",
        assets = ["assets/custom.css"],
    ),
    plugins = [bib],
    clean = true
)

# Version menu: stable (latest tag), main (main branch), dev (dev branch).
# Each deploydocs call only deploys when CI runs on its devbranch (or on a tag).
versions = ["stable" => "v^", "main" => "main", "dev" => "dev"]

deploydocs(
    repo = "github.com/awi-esc/ModelWeights.jl.git",
    devbranch = "main",
    devurl = "main",
    versions = versions
)

deploydocs(
    repo = "github.com/awi-esc/ModelWeights.jl.git",
    devbranch = "dev",
    devurl = "dev",
    versions = versions
)
