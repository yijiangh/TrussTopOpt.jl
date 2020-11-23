module TrussTopOpt

using Reexport

# Topology optimization problem definitions
include(joinpath("TrussTopOptProblems", "TrussTopOptProblems.jl"))

@reexport using .TrussTopOptProblems

end # module
