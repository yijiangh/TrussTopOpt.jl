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
    i = 1
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

    solver = FEASolver(Displacement, Direct, problem)
    # call solver to trigger assemble!
    solver()
    @show solver.u
    # scene, layout = draw_truss_problem(problem; u=solver.u)

    try
        global K, Kσ = buckling(problem, solver.globalinfo, solver.elementinfo)
    catch err
        println(RED_BG("ERROR: "))
        println(err)
        if isa(err, SingularException)
            println(RED_FG("Linear elasticity solve failed."))
            K_l = solver.globalinfo.K
            # TODO: use sparse eigen
            # λ, ϕ = eigs(K.data, nev = 2, which=:SM);
            # @show λ
            # @show ϕ

            # ? is this ordered by eigen values' magnitude?
            F = eigen(Array(K_l))
            @show F
            # @assert abs(F.values[1] - 0.0) < eps()

            # * draw eigen mode
            # scene, layout = draw_truss_problem(problem; u=F.vectors[:,1])
        end
    end
    println(GREEN_FG("Linear elasticity solved."))

    # Find the maximum eigenvalue of system Kσ x = 1/λ K x
    # https://julialinearalgebra.github.io/IterativeSolvers.jl/stable/eigenproblems/lobpcg/
    r = lobpcg(Kσ.data, K.data, true, 1)

    # # Minimum eigenvalue of the system K x = λ Kσ x
    @show λ = 1/r.λ[1]

    F = eigen(Array(Kσ), Array(K))
    # F = eigen(Array(K), Array(Kσ))
    @show 1 ./ F.values

    # if 0 < λ && λ < 1
        # 0 < λ_min < 1 means that the equilibrium equation with load 
        # λ*f is not solvable for some 0 < λ < 1 and the structure under load
        # f is not stable

        # s = 1.0
        # Ks = K .+ λ .* Kσ

        # eigen_vals, eigen_vecs = eigs(K, nev = 2, which=:SM);
        # @show eigen_vals
        # @show eigen_vecs
        # F = eigen(Array(K))
        # @show F

        # cholKs = cholesky(Symmetric(Ks))
        # u_b = cholKs \ solver.globalinfo.f
        # println(GREEN_FG("Buckling problem solved."))

        ## TODO plot analysis result with
        # scene, layout = draw_truss_problem(problem; u=u_b)
        # println("End")
    # end

    # TODO eigenmode to show global instability mode
    # TODO finite difference to verify the gradient
    # TODO verify the gradient for some analytical problems
    # TODO "manual" interior point loop, adjusting the c value every iter

# end # end test set