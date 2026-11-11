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

deploydocs(
    repo = "github.com/awi-esc/ModelWeights.jl.git",
    devbranch = "main"
)
