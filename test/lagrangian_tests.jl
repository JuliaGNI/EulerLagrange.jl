using EulerLagrange
using GeometricEquations
using LinearAlgebra
using ModelingToolkit
using Symbolics

@testset "Test Lagrangian" begin

    lvars = lagrangian_variables(2)

    @parameters t
    @variables (x(t))[1:2]
    @variables (v(t))[1:2]

    @test isequal(t, lvars[1])
    @test all([isequal(x[i], lvars[2][i]) for i in eachindex(x, lvars[2])])
    @test all([isequal(v[i], lvars[3][i]) for i in eachindex(v, lvars[3])])


    # Particle in square potential

    L = v ⋅ v / 2 - x ⋅ x / 2

    lag_sys = LagrangianSystem(t, x, v, L)


    # Lotka-Volterra System

    a₁ = -1.0
    a₂ = -1.0
    b₁ = 1.0
    b₂ = 2.0

    H = a₁ * x[1] + a₂ * x[2] + b₁ * log(x[1]) + b₂ * log(x[2])
    L = log(x[2]) / x[1] / 2 * v[1] - log(x[1]) / x[2] / 2 * v[2] - H

    lag_sys = LagrangianSystem(t, x, v, L)
    
end
