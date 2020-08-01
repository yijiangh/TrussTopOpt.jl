using Test
using BenchmarkTools

# using TopOpt.TopOptProblems
using TrussTopOpt

# E = 1.0
# ν = 0.3
# f = 1.0
ins_dir = joinpath(@__DIR__, "instances")
truss_file = joinpath(ins_dir, "tim.json")
load_supp_file = joinpath(ins_dir, "tim_load_support_case.json")

ndim, n, m, node_points, elements = parse_truss_json(truss_file)
loads, boundary = parse_support_load_json(load_supp_file)

tgrid = TrussTopOpt.TrussGrid(node_points, elements, loads, boundary)
