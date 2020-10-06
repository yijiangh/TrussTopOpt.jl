using Documenter, TrussTopOpt

makedocs(
    sitename = "TrussTopOpt.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "index.md",
        # "Quick start" => "QuickStart.md",
        # "TopOptProblems" => "TopOptProblems.md",
    ],
)

# if get(ENV, "CI", nothing) == "true"
#     deploydocs(
#         repo = "github.com/mohamed82008/TopOpt.jl.git",
#     )
# end

