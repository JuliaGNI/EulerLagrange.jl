using LinearAlgebra
using ModelingToolkit
using Symbolics

using Test
using EulerLagrange

@testset "Test Hamiltonian" begin

    t̄, q̄, p̄ = hamiltonian_variables(2)

    @parameters t
    q = @variables (q(t))[1:2]
    p = @variables (p(t))[1:2]

    @test isequal(t, t̄)
    @test all([isequal(q[i], q̄[i]) for i in eachindex(q,q̄)])
    @test all([isequal(p[i], p̄[i]) for i in eachindex(p,p̄)])


    H = p ⋅ p / 2 + q ⋅ q / 2

    ham_sys = HamiltonianSystem(t, q, p, H)
    
end
