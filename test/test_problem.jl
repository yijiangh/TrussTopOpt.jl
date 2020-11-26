using TopOpt
using TrussTopOpt
using Makie

ins_dir = joinpath(@__DIR__, "instances", "ground_meshes");

@testset "Tim problem to solve" for lc_ind in [0, 1]
    file_name = "tim.json"
    problem_file = joinpath(ins_dir, file_name)

    node_points, elements, Es, crosssecs, fixities, load_cases = parse_truss_json(problem_file);
    loads = load_cases[string(lc_ind)]

    problem = TrussProblem(Val{:Linear}, node_points, elements, loads, fixities, Es, crosssecs);

    ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)
    @test problem.E == Es

    penalty = TopOpt.PowerPenalty(1.0) # 3
    solver = FEASolver(Displacement, Direct, problem, xmin = xmin,
        penalty = penalty);

    # TopOpt.LogBarrier
    # linear_elasticity, du/dx
    # Compliance
    obj = Objective(TopOpt.Compliance(problem, solver, filterT = nothing,
        rmin = rmin, tracing = true, logarithm = false));

    constr = Constraint(TopOpt.Volume(problem, solver, filterT = nothing, rmin = rmin), V);

    mma_options = options = MMA.Options(maxiter = 3000,
        tol = MMA.Tolerances(kkttol = 0.001))
    convcriteria = MMA.KKTCriteria()
    optimizer = MMAOptimizer(obj, constr, MMA.MMA87(),
        ConjugateGradient(), options = mma_options,
        convcriteria = convcriteria);

    simp = SIMP(optimizer, penalty.p);

    # ? 1.0 might induce an infeasible solution, which gives the optimizer a hard time to escape 
    # from infeasible regions and return a result
    x0 = fill(0.5, length(solver.vars))
    result = simp(x0);

    scene, layout = draw_truss_problem(problem; crosssecs=result.topology)
    # display(scene)

end # end testset