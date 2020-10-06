using JuAFEM

##########################################

ndim = 2
cdim = 2
fdim = 1 # field dimension

# grid = generate_grid(Line, (2,), Vec{1}((-1.,)), Vec{1}((1.,)))
grid = generate_grid(Quadrilateral, (2, 2));
dh = DofHandler(grid)
push!(dh, :u, fdim) # Add a displacement field
close!(dh)

fieldnames(typeof(grid))
grid.nodes

# JuAFEM.default_interpolation(getcelltype(dh.grid))
geom_order = 1
quad_order = 2
refshape = JuAFEM.getrefshape(dh.field_interpolations[1])

# create a new reference shape RefCube
# overload reference_coordinates
# separate dims

interpolation_space = Lagrange{cdim, refshape, geom_order}()
quadrature_rule = QuadratureRule{cdim, refshape}(quad_order)
cv = CellScalarValues(quadrature_rule, interpolation_space)

cv.N
cv.dNdx
cv.dNdξ
cv.detJdV
cv.M
cv.dMdξ
cv.qr_weights

n_basefuncs = getnbasefunctions(cv)

celliterator = CellIterator(dh)
cells = collect(enumerate(celliterator))
reinit!(cv, cells[1][2])
cv

ci = cells[1][2]
ci.coords

getnquadpoints(cellvalues)
ci.coords[2]
cv.dMdξ[1, 1]
ci.coords[2] ⊗ cv.dMdξ[1, 1]

A = Tensor{1, 2}((1,2))
# At = Tensor{2, 1}((1),(2))
At = rand(Tensor{2, 1})
B = Tensor{1,1}((1,))
A ⊗ B
At ⊗ B

ci.coords[1] ⊗ cv.dMdξ[1, 1]

E = 1.
A = 1.
C = Tensor{1, 1}([E])
Kesize = ndofs_per_cell(dh, 1)
Ke_e = zeros(Float64, ndim, ndim)
fe = zeros(Float64, Kesize)
Ke_0 = Matrix{Float64}(undef, Kesize, Kesize)

for q_point in 1:getnquadpoints(cellvalues)
    @show dΩ = getdetJdV(cellvalues, q_point)
    for b in 1:n_basefuncs
        @show ∇ϕb = shape_gradient(cellvalues, q_point, b)
        @show ϕb = shape_value(cellvalues, q_point, b)
        for d2 in 1:ndim
            # fe = @set fe[(b-1)*dim + d2] += ϕb * body_force[d2] * dΩ
            for a in 1:n_basefuncs
                @show ∇ϕa = shape_gradient(cellvalues, q_point, a)
                Ke_e .= dotdot(∇ϕa, C, ∇ϕb) * dΩ
                # Ke_e .= ∇ϕa ⋅ C ⋅ ∇ϕb * dΩ * A
                for d1 in 1:ndim
                    #if dim*(b-1) + d2 >= dim*(a-1) + d1
                    Ke_0[ndim*(a-1) + d1, ndim*(b-1) + d2] += Ke_e[d1,d2]
                    #end
                end
            end
        end
    end
end

Ke_0

