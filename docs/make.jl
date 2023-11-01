using Documenter, SorptionModels

makedocs(;
    modules=[SorptionModels],
    authors="Will <william.joseph.box@gmail.com> and contributors",
    repo="https://github.com/Boxylmer/SorptionModels.jl/blob/{commit}{path}#{line}",
    sitename="SorptionModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Boxylmer.github.io/SorptionModels.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Models" => [
            "Empirical Models" => "models/empirical_models.md",
            "Fundamental Models" => "models/fundamental_models.md",
            "Transient Models" => "models/transient_models.md",
            "Dilation Models" => "models/dilation_models.md"
        ],
        "Model Analyses" => "model analyses.md",
        "Analysis Writers" => "writers.md",
        "Internals" => "internals.md"
    ],
)


deploydocs(;
    repo="github.com/Boxylmer/SorptionModels.jl.git",
    devbranch="master",
)