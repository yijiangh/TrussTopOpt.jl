module TrussTopOpt

using Reexport

# Topology optimization problem definitions
include(joinpath("TrussTopOptProblems", "TrussTopOptProblems.jl"))
@reexport using .TrussTopOptProblems

# Objective and constraint functions
include(joinpath("Functions", "Functions.jl"))
@reexport using .Functions

end # module
