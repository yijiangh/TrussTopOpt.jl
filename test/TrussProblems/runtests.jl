using Test
using BenchmarkTools

# using TopOpt.TopOptProblems
using TrussTopOpt
using JuAFEM
using Setfield

# E = 1.0
# ν = 0.3
# f = 1.0
ins_dir = joinpath(@__DIR__, "instances")
truss_file = joinpath(ins_dir, "tim.json")
load_supp_file = joinpath(ins_dir, "tim_load_support_case.json")

ndim, n, m, node_points, elements = parse_truss_json(truss_file)
loads, boundary = parse_support_load_json(load_supp_file)

tgrid = TrussGrid(node_points, elements, boundary)
# # @test JuAFEM.getnnodes(tgrid.grid) == n
# # @test JuAFEM.getncells(tgrid.grid) == m

fieldnames(typeof(tgrid))
getcoordinates(tgrid.grid.nodes[1])
nfaces(tgrid.grid.cells[1])
JuAFEM.nnodes(tgrid.grid.cells[1])

dim = typeof(tgrid.grid.cells[1]).parameters[1]

Line

@show JuAFEM.celltypes

# dh = DofHandler(tgrid.grid)
# dim = 2
# ip = Lagrange{dim,RefCube,1}()
# BCValues(ip, default_interpolation(getcelltype(dh.grid))
# push!(dh, :u, dim, ip)

# tproblem = TrussProblem(Val{:Linear}, node_points, elements, loads, boundary)

##########################################

dim = 1
line_grid = generate_grid(Line, (2,), Vec{1}((-2.,)), Vec{1}((2.,)))
dh = DofHandler(line_grid)
push!(dh, :u, dim) # Add a displacement field
close!(dh)

fieldnames(typeof(line_grid))
line_grid.nodes

# JuAFEM.default_interpolation(getcelltype(dh.grid))
geom_order = 1
quad_order = 2
refshape = JuAFEM.getrefshape(dh.field_interpolations[1])

interpolation_space = Lagrange{dim, refshape, geom_order}()
quadrature_rule = QuadratureRule{dim, refshape}(quad_order)
# interpolation_space = Lagrange{dim-1, refshape, geom_order}()
# quadrature_rule = QuadratureRule{dim-1, refshape}(quad_order)
cellvalues = CellScalarValues(quadrature_rule, interpolation_space)
facevalues = FaceScalarValues(QuadratureRule{dim-1, refshape}(quad_order), interpolation_space)

cellvalues.N
facevalues.N
n_basefuncs = getnbasefunctions(cellvalues)

celliterator = CellIterator(dh)
cells = collect(enumerate(celliterator))
reinit!(cellvalues, cells[1][2])

getnquadpoints(cellvalues)

E = 1.
A = 1.
C = Tensor{1, 1}([E])
Kesize = ndofs_per_cell(dh, 1)
Ke_e = zeros(Float64, dim, dim)
fe = zeros(Float64, Kesize)
Ke_0 = Matrix{Float64}(undef, Kesize, Kesize)

for q_point in 1:getnquadpoints(cellvalues)
    @show dΩ = getdetJdV(cellvalues, q_point)
    for b in 1:n_basefuncs
        @show ∇ϕb = shape_gradient(cellvalues, q_point, b)
        @show ϕb = shape_value(cellvalues, q_point, b)
        for d2 in 1:dim
            # fe = @set fe[(b-1)*dim + d2] += ϕb * body_force[d2] * dΩ
            for a in 1:n_basefuncs
                @show ∇ϕa = shape_gradient(cellvalues, q_point, a)
                # Ke_e .= dotdot(∇ϕa, C, ∇ϕb) * dΩ
                Ke_e .= ∇ϕa ⋅ C ⋅ ∇ϕb * dΩ * A
                for d1 in 1:dim
                    #if dim*(b-1) + d2 >= dim*(a-1) + d1
                    Ke_0[dim*(a-1) + d1, dim*(b-1) + d2] += Ke_e[d1,d2]
                    #end
                end
            end
        end
    end
end

Ke_0
Ke_e

##########################################
# beam element

dim = 1
side_len = 1.0
line_grid = generate_grid(Line, (2,), Vec{1}((-side_len,)), Vec{1}((side_len,)))
dh = DofHandler(line_grid)

push!(dh, :u, 4) # Add a displacement field
close!(dh)

# JuAFEM.default_interpolation(getcelltype(dh.grid))
geom_order = 3
quad_order = 2
refshape = JuAFEM.getrefshape(dh.field_interpolations[1])

interpolation_space = Lagrange{dim, refshape, geom_order}()
quadrature_rule = QuadratureRule{dim, refshape}(quad_order)
# interpolation_space = Lagrange{dim-1, refshape, geom_order}()
# quadrature_rule = QuadratureRule{dim-1, refshape}(quad_order)
cellvalues = CellScalarValues(quadrature_rule, interpolation_space)
facevalues = FaceScalarValues(QuadratureRule{dim-1, refshape}(quad_order), interpolation_space)

cellvalues.N
facevalues.N
n_basefuncs = getnbasefunctions(cellvalues)

celliterator = CellIterator(dh)
cells = collect(enumerate(celliterator))
reinit!(cellvalues, cells[1][2])

getnquadpoints(cellvalues)

E = 1.
A = 1.
C = Tensor{1, 1}([E])
Kesize = ndofs_per_cell(dh, 1)
Ke_e = zeros(Float64, dim, dim)
fe = zeros(Float64, Kesize)
Ke_0 = Matrix{Float64}(undef, Kesize, Kesize)

for q_point in 1:getnquadpoints(cellvalues)
    @show dΩ = getdetJdV(cellvalues, q_point)
    for b in 1:n_basefuncs
        @show ∇ϕb = shape_gradient(cellvalues, q_point, b)
        @show ϕb = shape_value(cellvalues, q_point, b)
        for d2 in 1:dim
            # fe = @set fe[(b-1)*dim + d2] += ϕb * body_force[d2] * dΩ
            for a in 1:n_basefuncs
                @show ∇ϕa = shape_gradient(cellvalues, q_point, a)
                # Ke_e .= dotdot(∇ϕa, C, ∇ϕb) * dΩ
                Ke_e .= ∇ϕa ⋅ C ⋅ ∇ϕb * dΩ * A
                for d1 in 1:dim
                    #if dim*(b-1) + d2 >= dim*(a-1) + d1
                    Ke_0[dim*(a-1) + d1, dim*(b-1) + d2] += Ke_e[d1,d2]
                    #end
                end
            end
        end
    end
end

Ke_0
Ke_e

##########################################

using TopOpt
using TimerOutputs
const to = TimerOutput()

E = 1.0 # Young’s modulus
v = 0.3 # Poisson’s ratio
f = 1.0 # downward force
problem = @timeit to "build problem" PointLoadCantilever(Val{:Linear}, (40, 20, 20), (1.0, 1.0, 1.0), E, v, f);
V = 0.3 # volume fraction
xmin = 0.001 # minimum density
rmin = 4.0 # density filter radius

penalty = TopOpt.PowerPenalty(3.0)
solver = @timeit to "build fea solve" FEASolver(Displacement, Direct, problem, xmin = xmin, penalty = penalty);

@timeit to "solve" solver();
solver.lhs
solver.rhs

to

# ke = E*A/Le*[ 1 -1; -1 1 ]
# shape function N = [(L-x)/L x/L]
# B = d/dx (N) = [-1/L 1/L]
# k = ∫_0^L B^T E B A dx = AE/L [1 -1; -1 1]
