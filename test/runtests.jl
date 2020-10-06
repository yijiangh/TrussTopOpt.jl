using Test, SafeTestsets

@safetestset "Problem Def Tests" begin include("test_problems.jl") end
# @safetestset "FEA Tests" begin include("test_fea.jl") end