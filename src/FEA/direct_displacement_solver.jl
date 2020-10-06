using TopOpt.FEA: DirectDisplacementSolver

function DirectDisplacementSolver(sp::TrussProblem{dim, T};
    xmin=T(1)/1000, 
    penalty=PowerPenalty{T}(1), 
    prev_penalty=copy(penalty),
    quad_order=default_quad_order(sp)) where {dim, T}

    elementinfo = ElementFEAInfo(sp, quad_order, Val{:Static})
    globalinfo = GlobalFEAInfo(sp)
    u = zeros(T, ndofs(sp.ch.dh))
    lhs = similar(u)
    rhs = similar(u)
    vars = fill(one(T), getncells(sp.ch.dh.grid) - sum(sp.black) - sum(sp.white))
    varind = sp.varind

    prev_penalty = setpenalty(prev_penalty, T(NaN))
    return DirectDisplacementSolver(sp, globalinfo, elementinfo, u, lhs, rhs, vars, penalty, prev_penalty, xmin)
end