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
    # i = 2
    file_name = problem_json[i]
    problem_file = joinpath(ins_dir, file_name)

    node_points, elements, Es, crosssecs, fixities, load_cases = parse_truss_json(problem_file);
    ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)
    loads = load_cases["0"]

    problem = TrussProblem(Val{:Linear}, node_points, elements, loads, fixities, Es, crosssecs);

    @test getdim(problem) == ndim
    @test JuAFEM.getncells(problem) == ncells
    @test problem.E == Es
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
        Ke_m = (A*Es[cellidx]/L)*Γ'*[1 -1; -1 1]*Γ
        Ke = elementinfo.Kes[cellidx]
        @test Ke_m ≈ Ke
    end

    solver = FEASolver(Displacement, Direct, problem)
    solver()

    # TODO plot analysis result with
    scene, layout = draw_truss_problem(problem)

    # we use kN for force and m for length
    # thus, pressure/modulus is in kN/m
    # the textbook uses a quite rough rounding scheme...
    # 0.3 mm error
    to_K_full = solver.globalinfo.K.data
    @assert norm(solver.u - u_solutions[i]) < 3e-4

end # end test set