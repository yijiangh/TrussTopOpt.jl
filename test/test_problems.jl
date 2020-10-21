using Test
using LinearAlgebra
using JuAFEM
using TopOpt
using TopOpt.TopOptProblems: boundingbox, nnodespercell, getgeomorder, getmetadata, getdh, getE, getdim

include("utils.jl")

# test cell volumes
using JuAFEM: cellid, getcoordinates, CellIterator
using TrussTopOpt.TrussTopOptProblems
using TrussTopOpt.TrussTopOptProblems: gettrussgrid, getJuaFEMgrid, getcrosssecs, default_quad_order, make_Kes_and_fes

# tim problem
# @testset "tim_problem" begin

ins_dir = joinpath(@__DIR__, "instances", "fea_examples");
truss_file = joinpath(ins_dir, "mgz_truss2.json");
load_supp_file = joinpath(ins_dir, "mgz_truss2_load_support.json");

ndim, nnodes, ncells, node_points, elements, E, crosssecs = parse_truss_json(truss_file);
loads, boundary = parse_support_load_json(load_supp_file);

problem = TrussProblem(Val{:Linear}, node_points, elements, loads, boundary, E, crosssecs);

@test getdim(problem) == ndim
@test JuAFEM.getncells(problem) == ncells
@test problem.E == E
@test problem.black == problem.white == falses(ncells)
@test problem.force == loads
@test problem.varind == 1:ncells
grid = problem.ch.dh.grid
@test length(grid.cells) == ncells

@test getgeomorder(problem) == 1
@test nnodespercell(problem) == 2

# ! specific for "tim.json"
# corners = [[0.0, 0.0], [8.5, 1.0]]
# for i in 1:2, j in 1:2
#     @test boundingbox(grid)[i][j] ≈ corners[i][j] atol = 1e-8
# end
# @test length(grid.boundary_matrix.nzval) == length(boundary) * 2

# metadata = getmetadata(problem)

quad_order=default_quad_order(problem)
elementinfo = ElementFEAInfo(problem, quad_order, Val{:Static});
for cell in CellIterator(getdh(problem))
    cellidx = cellid(cell)
    coords = getcoordinates(cell) # get the coordinates
    L = norm(coords[1] - coords[2])
    A = getcrosssecs(problem)[cellidx]
    @test elementinfo.cellvolumes[cellidx] ≈ L * A

    Γ = global2local_transf_matrix(coords...)
    book_Ke = (A*E[cellidx]/L)*Γ'*[1 -1; -1 1]*Γ
    Ke = elementinfo.Kes[cellidx]
    @test book_Ke ≈ Ke
end

solver = FEASolver(Displacement, Direct, problem)
solver()
solver.u

# end # end test set