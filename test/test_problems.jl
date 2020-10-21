using Test
using TopOpt
using TopOpt.TopOptProblems: boundingbox, nnodespercell, getgeomorder, getmetadata, getdh, getE
using TrussTopOpt.TrussTopOptProblems
import TrussTopOpt.TrussTopOptProblems: default_quad_order, make_Kes_and_fes
using JuAFEM
using LinearAlgebra: norm

E = 1.0
ν = 0.3
f = 1.0

# tim problem
# @testset "tim_problem" begin
# global E, ν, force

ins_dir = joinpath(@__DIR__, "instances");
# truss_file = joinpath(ins_dir, "tim.json");
# load_supp_file = joinpath(ins_dir, "tim_load_support_case.json");
truss_file = joinpath(ins_dir, "mgz_truss1.json");
load_supp_file = joinpath(ins_dir, "mgz_truss1_load_support.json");

ndim, nnodes, ncells, node_points, elements = parse_truss_json(truss_file);
loads, boundary = parse_support_load_json(load_supp_file);

problem = TrussProblem(Val{:Linear}, node_points, elements, loads, boundary);

# @test JuAFEM.getncells(problem) == ncells
# @test problem.E == E
# @test problem.ν == ν
# @test problem.black == problem.white == falses(ncells)
# @test problem.force == loads
# @test problem.varind == 1:ncells
# grid = problem.ch.dh.grid
# @test length(grid.cells) == ncells

# @test getgeomorder(problem) == 1
# @test nnodespercell(problem) == 2

# corners = [[0.0, 0.0], [8.5, 1.0]]
# for i in 1:2, j in 1:2
#     @test boundingbox(grid)[i][j] ≈ corners[i][j] atol = 1e-8
# end
# @test length(grid.boundary_matrix.nzval) == length(boundary) * 2

quad_order=default_quad_order(problem)
elementinfo = ElementFEAInfo(problem, quad_order, Val{:Static});
metadata = getmetadata(problem)

# test cell volumes
using JuAFEM: cellid, getcoordinates, CellIterator
using TrussTopOpt.TrussTopOptProblems: gettrussgrid, getJuaFEMgrid, getcrosssecs

function global2local_transf_matrix(end_vert_u, end_vert_v)
    @assert length(end_vert_u) == length(end_vert_v)
    @assert length(end_vert_u) == 2 || length(end_vert_u) == 3
    xdim = length(end_vert_u)
    L = norm(end_vert_u-end_vert_v)
    @assert L > 1e-6

    # by convention, the new x axis is along the element's direction
    # directional cosine of the new x axis in the global world frame
    c_x = (end_vert_v[1] - end_vert_u[1])/L
    c_y = (end_vert_v[2] - end_vert_u[2])/L
    R = zeros(2,xdim*2)

    if 3 == xdim
        error("Not Implemented")
        # c_z = (end_vert_v[2] - end_vert_u[2]) / L
        # # TODO rotaxis
        # if abs(abs(c_z) - 1.0) < eps()
        #     # the element is parallel to global z axis
        #     # cross product is not defined, in this case
        #     # it's just a rotation about the global z axis
        #     # in x-y plane
        #     R[1, 3] = -c_z
        #     R[2, 2] = 1
        #     R[3, 1] = c_z
        # else
        #     # local x_axis = element's vector
        #     new_x = [c_x, c_y, c_z]
        #     # local y axis = cross product with global z axis
        #     new_y = -cross(new_x, [0,0,1.0])
        #     new_y /= norm(new_y)
        #     new_z = cross(new_x, new_y)
        #     R[1, :] = new_x
        #     R[2, :] = new_y
        #     R[3, :] = new_z
        # end
    elseif 2 == xdim
        R = [c_x c_y 0 0; 0 0 c_x c_y]
    end
    return R
end

E = getE(problem)
for cell in CellIterator(getdh(problem))
    cellidx = cellid(cell)
    coords = getcoordinates(cell) # get the coordinates
    L = norm(coords[1] - coords[2])
    A = getcrosssecs(problem)[cellidx]
    @test elementinfo.cellvolumes[cellidx] ≈ L * A

    R = global2local_transf_matrix(coords...)
    book_Ke = (A*E/L)*R'*[1 -1; -1 1]*R
    Ke = elementinfo.Kes[cellidx]
    @test book_Ke ≈ Ke
end

# end # end test set