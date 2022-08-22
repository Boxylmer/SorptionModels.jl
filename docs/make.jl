using Documenter, SorptionModels

makedocs(
    sitename="SorptionModels.jl",
    modules = [SorptionModels],
    format = Documenter.HTML(prettyurls = false),
)

# deploydocs(
#     repo = "github.com/Boxylmer/SorptionModels.jl.git",
# )

# deploydocs(
#     target = "build",
#     repo = "github.com/Boxylmer/SorptionModels.jl.git",
#     deps   = nothing,
#     make   = nothing,
#     push_preview = true
# )