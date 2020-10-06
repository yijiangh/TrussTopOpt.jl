using Test

###########################################
# check init cellvalues for truss element

using Tensors
using JuAFEM
using TrussTopOpt
const TTO = TrussTopOpt

T = Float64
refshape = RefCube
# nodal embedding dimension
xdim = 2

# cell dimension
# i.e. ξ dimension
cdim = 1
# geometric interpolation order
geom_order = 1

# quadrature integration order
quad_order = 2

func_interpol = Lagrange{cdim, refshape, geom_order}()
quad_rule = QuadratureRule{cdim, refshape}(quad_order)
# geom_interpol = Lagrange{cdim, refshape, geom_order}()
# geom_quad_rule = QuadratureRule{cdim, refshape}(quad_order)

# cv = CellScalarValues(quad_rule, func_interpol)
cv = CellScalarValues(T, quad_rule, func_interpol; xdim=xdim)

x = [Vec{ndim}((1.,1)), Vec{ndim}((2.,2.))]
cv.dNdξ
cv.dNdx

A = 1.0

TrussTopOptProblems.truss_reinit!(cv, x, A)
cv.dNdξ
cv.dNdx

# L = norm(x[1]-x[2])
# 1/L
# Line

##############################
# ideal output

# ke = E*A/Le*[ 1 -1; -1 1 ]
# shape function N = [(L-x)/L x/L]
# B = d/dx (N) = [-1/L 1/L]
# k = ∫_0^L B^T E B A dx = AE/L [1 -1; -1 1]

# end