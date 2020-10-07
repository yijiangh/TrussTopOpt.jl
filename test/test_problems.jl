using Test
using TopOpt
using TopOpt.TopOptProblems: boundingbox, nnodespercell, getgeomorder
using TrussTopOpt.TrussTopOptProblems
import TrussTopOpt.TrussTopOptProblems: default_quad_order, make_Kes_and_fes
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
@test nnodespercell(problem) == 2

corners = [[0.0, 0.0], [8.5, 1.0]]
for i in 1:2, j in 1:2
    @test boundingbox(grid)[i][j] ≈ corners[i][j] atol = 1e-8
end
@test length(grid.boundary_matrix.nzval) == length(boundary) * 2

quad_order=default_quad_order(problem)

fnode = Tuple(getnodeset(problem.truss_grid.grid, "load"))[1]
problem.force

elementinfo = ElementFEAInfo(problem, quad_order, Val{:Static})

Kes, weights, cellvalues = make_Kes_and_fes(problem, quad_order, Val{:Static})

element_Kes = convert(
    Vector{<:ElementMatrix},
    Kes;
    bc_dofs = problem.ch.prescribed_dofs,
    dof_cells = problem.metadata.dof_cells,
)


# end # end test set