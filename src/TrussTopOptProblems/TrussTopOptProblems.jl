module TrussTopOptProblems

using JuAFEM, StaticArrays, LinearAlgebra
using SparseArrays
using TopOpt
using TopOpt.Utilities
using Setfield
# using ..TopOpt: PENALTY_BEFORE_INTERPOLATION
# using ..Utilities: @forward_property

using VTKDataTypes
#using Makie
#using GeometryTypes

import JuAFEM: assemble!

# include("utils.jl")
include("grids.jl")
# include("metadata.jl")
include("problem_types.jl")
include("matrices_and_vectors.jl")
include("elementinfo.jl")
include(joinpath("TrussIO", "TrussIO.jl"))
using .TrussIO
include(joinpath("TrussPlotting", "TrussPlotting.jl"))
using .TrussPlotting

export TrussGrid, TrussProblem
export parse_truss_json
export draw_truss_problem!, draw_truss_problem

end # module
