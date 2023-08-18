using EulerLagrange
using EulerLagrange: equations
using GeometricEquations
using LinearAlgebra
using ModelingToolkit
using Test


@testset "Test Hamiltonian" begin

    hvars = hamiltonian_variables(2)

    @parameters t
    @variables (q(t))[1:2]
    @variables (p(t))[1:2]

    @test isequal(t, hvars[1])
    @test all([isequal(q[i], hvars[2][i]) for i in eachindex(q, hvars[2])])
    @test all([isequal(p[i], hvars[3][i]) for i in eachindex(p, hvars[3])])


    H(t, q, p) = p ⋅ p / 2 + q ⋅ q / 2

    ham_sys = HamiltonianSystem(H(t, q, p), t, q, p)

    ṽ(v, t, q, p) = v .= p
    f̃(f, t, q, p) = f .= -q
    ż̃(ż, t, q, p) = ż .= [p..., -q...]

    t₀, q₀, p₀ = (0.0, [1.0, 1.0], [0.5, 2.0])
    z₀ = [q₀..., p₀...]

    v₁, v₂ = zero(q₀), zero(q₀)
    f₁, f₂ = zero(p₀), zero(p₀)
    ż₁, ż₂ = zero(z₀), zero(z₀)
        
    eqs = equations(ham_sys)

    eqs.v(v₁, t₀, q₀, p₀)
    eqs.f(f₁, t₀, q₀, p₀)
    eqs.ż(ż₁, t₀, q₀, p₀)

    ṽ(v₂, t₀, q₀, p₀)
    f̃(f₂, t₀, q₀, p₀)
    ż̃(ż₂, t₀, q₀, p₀)
        

    @test eqs.H(t₀, q₀, p₀) == H(t₀, q₀, p₀)
    @test v₁ == v₂
    @test f₁ == f₂
    @test ż₁ == ż₂

end
