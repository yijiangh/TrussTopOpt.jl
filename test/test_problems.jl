using Test

using TopOpt
using TopOpt.TopOptProblems: boundingbox, nnodespercell
using TrussTopOpt.TrussTopOptProblems
# import TrussTopOpt.TrussTopOptProblems: nnodespercell
using JuAFEM

E = 1.0
ν = 0.3
f = 1.0

# tim problem
# @testset "tim_problem" begin
# global E, ν, force

ins_dir = joinpath(@__DIR__, "instances");
truss_file = joinpath(ins_dir, "tim.json");
load_supp_file = joinpath(ins_dir, "tim_load_support_case.json");

ndim, nnodes, ncells, node_points, elements = parse_truss_json(truss_file);
loads, boundary = parse_support_load_json(load_supp_file);

problem = TrussProblem(Val{:Linear}, node_points, elements, loads, boundary);

@test JuAFEM.getncells(problem) == ncells
@test problem.E == E
@test problem.ν == ν
@test problem.black == problem.white == falses(ncells)
@test problem.force == loads
@test problem.varind == 1:ncells
grid = problem.ch.dh.grid
@test length(grid.cells) == ncells

@test getgeomorder(problem) == 1

nnodespercell(problem)
@test nnodespercell(problem) == 2

corners = [[0.0, 0.0], [8.5, 1.0]]
for i in 1:2, j in 1:2
    @test boundingbox(grid)[i][j] ≈ corners[i][j] atol = 1e-8
end
@test length(grid.boundary_matrix.nzval) == length(boundary) * 2


# end # end test set