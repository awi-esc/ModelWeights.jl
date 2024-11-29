using Documenter, SimilarityWeights

makedocs(
    sitename = "SimilarityWeights.jl",
    modules = [SimilarityWeights],
    format   = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "index.md",
        "requirements.md",
        "getting-started.md",
        "weights.md",
        "examples/climwip.md",
        "examples/lgm.md",
        "references.md"
    ]
)