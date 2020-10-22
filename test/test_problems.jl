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
truss_file = joinpath(ins_dir, "mgz_truss3.json");
load_supp_file = joinpath(ins_dir, "mgz_truss3_load_support.json");

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
book_Ks = []
push!(book_Ks, ([0 -0.5 -0.5 0.5; 0 0 0.5 -0.5; 0 0 0 -0.5; 0 0 0 0], [0.5, 0.5, 0.5, 0.5], 707.11*1e3))
push!(book_Ks, ([0 0.433 -0.75 -0.433; 0 0 -0.433 -0.25; 0 0 0 0.433; 0 0 0 0], [0.75, 0.25, 0.75, 0.25], 375.00*1e3))
push!(book_Ks, ([0 0 -1 0; 0 0 0 0; 0 0 0 0; 0 0 0 0], [1.0, 0, 1.0, 0.0], 329.43*1e3))

quad_order=default_quad_order(problem)
elementinfo = ElementFEAInfo(problem, quad_order, Val{:Static});
for cell in CellIterator(getdh(problem))
    cellidx = cellid(cell)
    println("---\n Cell #$cellidx")

    coords = getcoordinates(cell)
    L = norm(coords[1] - coords[2])
    A = getcrosssecs(problem)[cellidx]
    @test elementinfo.cellvolumes[cellidx] ≈ L * A

    Γ = global2local_transf_matrix(coords...)
    book_Ke = (A*E[cellidx]/L)*Γ'*[1 -1; -1 1]*Γ
    Ke = elementinfo.Kes[cellidx]
    @test book_Ke ≈ Ke

    # @show (book_Ke .- matrix_from_upper(book_Ks[cellidx]...)) ./ book_Ke
end

# we use kN for force and m for length
# thus, pressure/modulus is in kN/m
solver = FEASolver(Displacement, Direct, problem)
solver()

# truss2
# dof_from_element = [1 2 3 4; 1 2 5 6]
# dof_stat = [0,0,1,1,1,1]

# truss3
dof_from_element = [1 2 3 4; 1 2 5 6; 3 4 5 6]
dof_stat = [0,0,0,1,1,1]

n_fixed_dof = sum(dof_stat)
n_free_dof = length(dof_stat)-n_fixed_dof
# force vector free dof portion
@show P_f = solver.globalinfo.f[1:n_free_dof]

# * TO solver result
to_K_full = solver.globalinfo.K.data
@show to_K_ff = Array(to_K_full[1:n_free_dof,1:n_free_dof])
@show to_u = solver.u[1:n_free_dof]
println("---")

# * directly from book
# book_K = [0 -0.1124 -1.7612
#           0 0 0.2103
#           0 0 0]
# book_K += book_K'
# book_K += Diagonal([2.522, 1.147, 8.757])
# ! 1 Mpa = 1000 Kn/m
# ! 200,000 Mpa = 2e8 Kn/m
# ! the book_k above is in `mm`, * 1e-3 to make it in unit `m`
# book_K *= 2e8 * 1e-3 
# book_u = book_K \ [200, 600, -800]
# @show book_K
# @show book_u
println("---")

# * manual construction
build_K_full = assemble_global_stiffness_matrix(elementinfo.Kes, nnodes, dof_from_element)
Perm = compute_permutation(dof_stat)
build_K = Perm*build_K_full*Perm'
@show build_K_ff = Array(build_K[1:n_free_dof,1:n_free_dof])
@show build_u = build_K_ff \ P_f
@assert build_K_ff*build_u ≈ P_f
println("---")

# @assert book_u ≈ solver.u[1:3]
# @assert book_K ≈ solver.globalinfo.K.data[1:3, 1:3]

# end # end test set