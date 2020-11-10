using Test
using LinearAlgebra
using JuAFEM
using JuAFEM: cellid, getcoordinates, CellIterator

using TopOpt
using TopOpt.TopOptProblems: boundingbox, nnodespercell, getgeomorder, getmetadata, getdh, getE, getdim
using TrussTopOpt.TrussTopOptProblems
using TrussTopOpt.TrussTopOptProblems: getcrosssecs, default_quad_order

include("utils.jl")

problem_json = ["mgz_truss1.json", "mgz_truss2.json", "mgz_truss3.json"]
u_solutions = [
    1e-3 * vcat([2.41, 0.72], zeros(2*2)),
    1e-3 * vcat([0.1783, 2.7222, -0.4863], zeros(4*3)),
    1e-3 * vcat([0.871, 1.244], [0,0], [-0.193, 0]),
]
ins_dir = joinpath(@__DIR__, "instances", "fea_examples");

@testset "Truss problem solve - $(problem_json[i])" for i in 1:length(problem_json)
    file_name = problem_json[i]
    truss_file = joinpath(ins_dir, file_name)
    load_supp_file = joinpath(ins_dir, split(file_name, ".json")[1]*"_load_support.json");

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

    quad_order=default_quad_order(problem)
    elementinfo = ElementFEAInfo(problem, quad_order, Val{:Static});
    for cell in CellIterator(getdh(problem))
        cellidx = cellid(cell)
        coords = getcoordinates(cell)
        L = norm(coords[1] - coords[2])
        A = getcrosssecs(problem)[cellidx]
        @test elementinfo.cellvolumes[cellidx] ≈ L * A

        Γ = global2local_transf_matrix(coords...)
        Ke_m = (A*E[cellidx]/L)*Γ'*[1 -1; -1 1]*Γ
        Ke = elementinfo.Kes[cellidx]
        @test Ke_m ≈ Ke
    end

    solver = FEASolver(Displacement, Direct, problem)
    solver()


    # we use kN for force and m for length
    # thus, pressure/modulus is in kN/m
    # the textbook uses a quite rough rounding scheme...
    # 0.3 mm error
    to_K_full = solver.globalinfo.K.data
    @assert norm(solver.u - u_solutions[i]) < 3e-4

end # end test set

    # n_fixed_dof = sum(dof_stat)
    # n_free_dof = length(dof_stat)-n_fixed_dof

    # # * manual construction
    # build_K_full = assemble_global_stiffness_matrix(elementinfo.Kes, nnodes, dof_from_element)
    # Perm = compute_permutation(dof_stat)
    # build_K = Perm*build_K_full*Perm'
    # @show build_K_ff = Array(build_K[1:n_free_dof,1:n_free_dof])
    # @show build_u = build_K_ff \ P_f
    # @assert build_K_ff*build_u ≈ P_f
    # println("---")

    # @assert book_u ≈ solver.u[1:3]
    # @assert book_K ≈ solver.globalinfo.K.data[1:3, 1:3]

    # ! specific for "tim.json"
    # corners = [[0.0, 0.0], [8.5, 1.0]]
    # for i in 1:2, j in 1:2
    #     @test boundingbox(grid)[i][j] ≈ corners[i][j] atol = 1e-8
    # end
    # @test length(grid.boundary_matrix.nzval) == length(boundary) * 2

    # metadata = getmetadata(problem)
    # truss2
    # dof_from_element = [1 2 3 4; 1 2 5 6]
    # dof_stat = [0,0,1,1,1,1]

    # truss3
    # dof_from_element = [1 2 3 4; 1 2 5 6; 3 4 5 6]
    # dof_stat = [0,0,0,1,1,1]

