module TrussTopOpt

using Reexport

# Topopology optimization problem definitions
include(joinpath("TrussTopOptProblems", "TrussTopOptProblems.jl"))

@reexport using .TrussTopOptProblems

end # module
