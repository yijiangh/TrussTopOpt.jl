using Test
using LinearAlgebra
using JuAFEM
using JuAFEM: cellid, getcoordinates, CellIterator

using TopOpt
using TopOpt.TopOptProblems: boundingbox, nnodespercell, getgeomorder, getmetadata, getdh, getE, getdim
using TrussTopOpt.TrussTopOptProblems
using TrussTopOpt.TrussTopOptProblems: buckling
using TrussTopOpt.TrussTopOptProblems: getA, default_quad_order
using IterativeSolvers
using Arpack
using Crayons.Box

problem_json = ["buckling_2d_global_instab.json", "buckling_2d_nodal_instab.json"]
ins_dir = joinpath(@__DIR__, "instances", "fea_examples");

# @testset "Buckling problem solve - $(problem_json[i])" for i in 1:length(problem_json)
    i = 2
    file_name = problem_json[i]
    problem_file = joinpath(ins_dir, file_name)

    node_points, elements, mats, crosssecs, fixities, load_cases = parse_truss_json(problem_file);
    ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)
    loads = load_cases["0"]

    problem = TrussProblem(Val{:Linear}, node_points, elements, loads, fixities, mats, crosssecs);

    @test getdim(problem) == ndim
    @test JuAFEM.getncells(problem) == ncells
    @test getE(problem) == [m.E for m in mats]
    @test problem.black == problem.white == falses(ncells)
    @test problem.force == loads
    @test problem.varind == 1:ncells
    grid = problem.ch.dh.grid
    @test length(grid.cells) == ncells

    @test getgeomorder(problem) == 1
    @test nnodespercell(problem) == 2

    # quad_order=default_quad_order(problem)
    # elementinfo = ElementFEAInfo(problem, quad_order, Val{:Static});
    solver = FEASolver(Displacement, Direct, problem)
    ## call solver to trigger assemble!
    try
        solver()
        @show solver.u
        # scene, layout = draw_truss_problem(problem; u=solver.u)

        K, Kσ = buckling(problem, solver.globalinfo, solver.elementinfo)
        println(GREEN_FG("Linear elasticity solved."))

        # Find the maximum eigenvalue of system Kσ x = 1/λ K x
        # https://julialinearalgebra.github.io/IterativeSolvers.jl/stable/eigenproblems/lobpcg/
        r = lobpcg(Kσ.data, K.data, true, 1)

        # Minimum eigenvalue of the system K x = λ Kσ x
        @show λ = 1/r.λ[1]

        s = 1.0
        Ks = K .+ s .* Kσ
        cholKs = cholesky(Symmetric(Ks))
        u_b = cholKs \ s .* solver.globalinfo.f
        println(GREEN_FG("Buckling problem solved."))

        ## TODO plot analysis result with
        # scene, layout = draw_truss_problem(problem; u=solver.u)
        scene, layout = draw_truss_problem(problem; u=u_b)
    catch err
        println(RED_FG("Linear elasticity failed."))
        println(err)
        K = solver.globalinfo.K
        # TODO: use sparse eigen
        # λ, ϕ = eigs(K.data, nev = 2, which=:SM);
        # @show λ
        # @show ϕ

        # ? is this ordered by eigen values' magnitude?
        F = eigen(Array(K))
        # @show F
        @assert abs(F.values[1] - 0.0) < eps()

        # * draw eigen mode
        scene, layout = draw_truss_problem(problem; u=F.vectors[:,1])
    end

    # TODO eigenmode to show global instability mode
    # TODO finite difference to verify the gradientj
    # TODO verify the gradient for some analytical problems
    # TODO "manual" interior point loop, adjusting the c value every iter

# end # end test set