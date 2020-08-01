using TopOpt.TopOptProblems
using Test
using BenchmarkTools

# E = 1.0
# ν = 0.3
# f = 1.0
HERE = println(@__DIR__)

# @btime problem = PointLoadCantilever(Val{:Linear}, (40, 20, 20), (1.0, 1.0, 1.0), E, ν, f)
# @test true == true