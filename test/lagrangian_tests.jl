using EulerLagrange
using EulerLagrange: equations
using GeometricEquations
using LinearAlgebra
using ModelingToolkit
using Test


@testset "Test Lagrangian" begin

    lvars = lagrangian_variables(2)

    @parameters t
    @variables (x(t))[1:2]
    @variables (v(t))[1:2]

    @test isequal(t, lvars[1])
    @test all([isequal(x[i], lvars[2][i]) for i in eachindex(x, lvars[2])])
    @test all([isequal(v[i], lvars[3][i]) for i in eachindex(v, lvars[3])])


    # Particle in square potential

    L(t,x,v) = v ⋅ v / 2 - x ⋅ x / 2

    lag_sys = LagrangianSystem(L(t, x, v), t, x, v)

    p̃(p, t, q, q̇) = p .= q̇
    f̃(f, t, q, q̇) = f .= -q
 
    t₀, q₀, v₀ = (0.0, [1.0, 1.0], [0.5, 2.0])
    p₀ = zero(v₀)

    p₁, p₂ = zero(p₀), zero(p₀)
    f₁, f₂ = zero(p₀), zero(p₀)
        
    eqs = equations(lag_sys)

    eqs.ϑ(p₁, t₀, q₀, v₀)
    eqs.f(f₁, t₀, q₀, v₀)

    p̃(p₂, t₀, q₀, v₀)
    f̃(f₂, t₀, q₀, v₀)

    @test eqs.L(t₀, q₀, v₀) == L(t₀, q₀, v₀)
    @test p₁ == p₂
    @test f₁ == f₂


    ntime = 1000
    tstep = 0.01
    tspan = (0.0, ntime * tstep)
    ics   = (q = q₀, p = p₀, λ = zero(q₀))

    lode = LODE(lag_sys)

    lprob1 = LODEProblem(lag_sys, tspan, tstep, ics)
    lprob2 = LODEProblem(lag_sys, tspan, tstep, q₀, p₀)

    @test lode == equation(lprob1) == equation(lprob2)


    # Lotka-Volterra System

    a₁ = -1.0
    a₂ = -1.0
    b₁ = 1.0
    b₂ = 2.0

    H_LV(t,x,v) = a₁ * x[1] + a₂ * x[2] + b₁ * log(x[1]) + b₂ * log(x[2])
    L_LV(t,x,v) = log(x[2]) / x[1] / 2 * v[1] - log(x[1]) / x[2] / 2 * v[2] - H_LV(t,x,v)

    lag_sys = LagrangianSystem(L_LV(t,x,v), t, x, v)

end
