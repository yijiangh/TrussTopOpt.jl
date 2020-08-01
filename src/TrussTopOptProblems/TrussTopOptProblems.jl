module TrussTopOptProblems

using JuAFEM, StaticArrays, LinearAlgebra
using TopOpt
using SparseArrays
# using ..TopOpt.Utilities
# using ..TopOpt: PENALTY_BEFORE_INTERPOLATION
# using ..Utilities: @forward_property

using VTKDataTypes
#using Makie
#using GeometryTypes

import JuAFEM: assemble!

# abstract type AbstractTopOptProblem end

# include("utils.jl")
include("grids_truss.jl")
# include("metadata.jl")
# include("problem_types.jl")
# include("multiload.jl")
# ! include("matrices_and_vectors.jl")
# include("assemble.jl")
include(joinpath("IO", "IO.jl"))
using .IO
# include("makie.jl")
include("makie_truss.jl")

# export PointLoadCantilever, HalfMBB, LBeam, TieBeam,
#     InpStiffness, StiffnessTopOptProblem, AbstractTopOptProblem,
#     GlobalFEAInfo, ElementFEAInfo, YoungsModulus, assemble, assemble_f!,
#     RaggedArray, ElementMatrix, rawmatrix, bcmatrix, save_mesh, RandomMagnitude, MultiLoad
export TrussGrid
export parse_truss_json, parse_support_load_json

end # module
