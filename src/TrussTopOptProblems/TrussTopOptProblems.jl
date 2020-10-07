module TrussTopOptProblems

using JuAFEM, StaticArrays, LinearAlgebra
using SparseArrays
using TopOpt
using TopOpt.Utilities
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
include(joinpath("IO", "IO.jl"))
using .IO
include("makie.jl")

# export PointLoadCantilever, HalfMBB, LBeam, TieBeam,
#     InpStiffness, StiffnessTopOptProblem, AbstractTopOptProblem,
#     GlobalFEAInfo, ElementFEAInfo, YoungsModulus, assemble, assemble_f!,
#     RaggedArray, ElementMatrix, rawmatrix, bcmatrix, save_mesh, RandomMagnitude, MultiLoad
export TrussGrid, TrussProblem # _LinearTrussGrid
export parse_truss_json, parse_support_load_json

end # module
