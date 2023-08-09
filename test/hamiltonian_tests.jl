using EulerLagrange
using LinearAlgebra
using ModelingToolkit
using Symbolics

@testset "Test Hamiltonian" begin

    hvars = hamiltonian_variables(2)

    @parameters t
    @variables (q(t))[1:2]
    @variables (p(t))[1:2]

    @test isequal(t, hvars[1])
    @test all([isequal(q[i], hvars[2][i]) for i in eachindex(q, hvars[2])])
    @test all([isequal(p[i], hvars[3][i]) for i in eachindex(p, hvars[3])])


    H = p ⋅ p / 2 + q ⋅ q / 2

    # ham_sys = HamiltonianSystem(t, q, p, H)
    
end
