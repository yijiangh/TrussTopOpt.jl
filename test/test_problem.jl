using TopOpt
using TrussTopOpt
using Makie

ins_dir = joinpath(@__DIR__, "instances", "ground_meshes");
# ins_dir = joinpath(@__DIR__, "instances", "fea_examples");

# @testset "Tim problem to solve"
    file_name = "tim.json"
    # file_name = "mgz_truss1.json"
    problem_file = joinpath(ins_dir, file_name)

    ndim, nnodes, ncells, node_points, elements, E, crosssecs = parse_truss_json(truss_file);
    loads, boundary = parse_support_load_json(load_supp_file);

    problem = TrussProblem(Val{:Linear}, node_points, elements, loads, boundary, E, crosssecs);

    # scene, layout = layoutscene(resolution = (1200, 900))
    # draw_truss_problem!(scene, layout, problem)
    # display(scene)

    # E = 1.0 # Young’s modulus
    # v = 0.3 # Poisson’s ratio
    # f = 1.0; # downward force

    # nels = (12, 6, 6) # change to (40, 20, 20) for a more high-res result
    # problem = PointLoadCantilever(Val{:Linear}, nels, (1.0, 1.0, 1.0), E, v, f);

    V = 0.3 # volume fraction
    xmin = 0.001 # minimum density
    rmin = 4.0; # density filter radius

    penalty = TopOpt.PowerPenalty(1.0) # 3
    solver = FEASolver(Displacement, Direct, problem, xmin = xmin,
        penalty = penalty);

    # TODO plot analysis result

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

    # ? 1.0 might induce an infeasible solution, which gives the optimizer a hard time to escape from infeasible regions and return a result
    x0 = fill(0.5, length(solver.vars))
    result = simp(x0);

    # TODO: use draw_truss!
    # result_mesh = GeometryBasics.Mesh(problem, result.topology);
    # mesh(result_mesh);

# end # end testset