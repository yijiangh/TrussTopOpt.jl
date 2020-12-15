using LinearAlgebra: logdet

abstract type AbstractMatrixFunction{T} <: AbstractFunction{T} end
abstract type AbstractSymmetricMatrixFunction{T} <: AbstractMatrixFunction{T} end

"""
LogDet(f(x)=A)
"""
@params mutable struct LogDet{T} <: AbstractMatrixFunction{T}
    f::AbstractSymmetricMatrixFunction{T} # f(x) -> A <: Symmetric
    grad::AbstractVector{T}
    fevals::Int
    maxfevals::Int
end

# function Base.logdet(f::AbstractMatrixFunction{T}) where {T}
#     return LogDet(f, similar(f.grad, T), 0, 10^8)
# end

TopOpt.dim(l::LogDet) = dim(l.f)
@inline function Base.getproperty(vf::LogDet, f::Symbol)
    f === :f && return getfield(vf, :f)
    f === :grad && return getfield(vf, :grad)
    f === :fevals && return getfield(vf, :fevals)
    f === :maxfevals && return getfield(vf, :maxfevals)
    return getproperty(getfield(vf, :f), f)
end
@inline function Base.setproperty!(vf::LogDet, f::Symbol, v)
    f === :f && return setfield!(vf, :f, v)
    f === :grad && return setfield!(vf, :grad, v)
    f === :fevals && return setfield!(vf, :fevals, v)
    f === :maxfevals && return setfield!(vf, :maxfevals, v)
    return setproperty!(getfield(vf, :f), f, v)
end

function (v::LogDet{T})(x, grad = v.grad) where {T}
    v.fevals += 1
    A = v.f(x) # .+ sqrt(eps(T))
    # v.f.grad is a Matrix too
    grad .= dot(inv(A), v.f.grad)
    if grad !== v.grad
        v.grad .= grad
    end
    return logdet(cholesky(A))
end