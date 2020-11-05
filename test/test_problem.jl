using TopOpt
using TrussTopOpt
using Makie

ins_dir = joinpath(@__DIR__, "instances", "ground_meshes");
# ins_dir = joinpath(@__DIR__, "instances", "fea_examples");

# @testset "Tim problem to solve"
    file_name = "tim.json"
    # file_name = "mgz_truss1.json"
    truss_file = joinpath(ins_dir, file_name)
    load_supp_file = joinpath(ins_dir, split(file_name, ".json")[1]*"_load_support.json");

    ndim, nnodes, ncells, node_points, elements, E, crosssecs = parse_truss_json(truss_file);
    loads, boundary = parse_support_load_json(load_supp_file);

    problem = TrussProblem(Val{:Linear}, node_points, elements, loads, boundary, E, crosssecs);

    scene = Scene()
    draw_truss_problem!(scene, problem)

    # V = 0.3 # volume fraction
    # xmin = 0.001 # minimum density
    # rmin = 4.0; # density filter radius

    # penalty = TopOpt.PowerPenalty(3.0)
    # solver = FEASolver(Displacement, Direct, problem, xmin = xmin,
    #     penalty = penalty);

    # obj = Objective(TopOpt.Compliance(problem, solver, filterT = DensityFilter,
    #     rmin = rmin, tracing = true, logarithm = false));

    # constr = Constraint(TopOpt.Volume(problem, solver, filterT = DensityFilter, rmin = rmin), V);

    # mma_options = options = MMA.Options(maxiter = 3000,
    #     tol = MMA.Tolerances(kkttol = 0.001))
    # convcriteria = MMA.KKTCriteria()
    # optimizer = MMAOptimizer(obj, constr, MMA.MMA87(),
    #     ConjugateGradient(), options = mma_options,
    #     convcriteria = convcriteria);

    # simp = SIMP(optimizer, penalty.p);

    # x0 = fill(1.0, length(solver.vars))
    # result = simp(x0);

    # TODO: use draw_truss!
    # result_mesh = GeometryBasics.Mesh(problem, result.topology);
    # mesh(result_mesh);
# end