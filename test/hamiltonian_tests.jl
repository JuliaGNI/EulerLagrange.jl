using EulerLagrange
using LinearAlgebra
using ModelingToolkit
using Symbolics

using Test
using EulerLagrange

@testset "Test Hamiltonian" begin

    hvars = hamiltonian_variables(2)

    @parameters t
    @variables (q(t))[1:2]
    @variables (p(t))[1:2]

    @test isequal(t, hvars[1])
    @test all([isequal(q[i], hvars[2][i]) for i in eachindex(q, hvars[2])])
    @test all([isequal(p[i], hvars[3][i]) for i in eachindex(p, hvars[3])])

    H = p ⋅ p / 2 + q ⋅ q / 2

    ham_sys = HamiltonianSystem(t, q, p, H)


    H̃(t,q,p) = (p ⋅ p)/2 + (q ⋅ q)/2
    f̃(t,q,p) = -q
    ṽ(t,q,p) = p
    ż̃(t,q,p) = [-p...,q...]

    t̃, q̃, p̃ = (1, [0.5, 0.4], [0.8, 0.9])
        
    eqs = ham_sys.eqs

    @test eqs.H(t̃, q̃, p̃) == H̃(t̃, q̃, p̃)
    @test eqs.f(t̃, q̃, p̃) == f̃(t̃, q̃, p̃)
    @test eqs.v(t̃, q̃, p̃) == ṽ(t̃, q̃, p̃)
    @test eqs.ż(t̃, q̃, p̃) == ż̃(t̃, q̃, p̃)

end