using LinearAlgebra
using ModelingToolkit
using Symbolics

@testset "Test Lagrangian" begin

    t̄, x̄, v̄ = lagrangian_variables(2)

    @parameters t
    x = @variables (x(t))[1:2]
    v = @variables (v(t))[1:2]

    @test isequal(t, t̄)
    @test all([isequal(x[i], x̄[i]) for i in eachindex(x,x̄)])
    @test all([isequal(v[i], v̄[i]) for i in eachindex(v,v̄)])


    L = v ⋅ v / 2 - x ⋅ x / 2

    lag_sys = LagrangianSystem(t, x, v, L)
    
end
