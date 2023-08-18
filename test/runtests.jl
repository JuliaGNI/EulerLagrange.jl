using EulerLagrange
using SafeTestsets

@safetestset HamiltonianTests = "$(rpad("Hamiltonian Systems",80))" begin include("hamiltonian_tests.jl") end
@safetestset LagrangianTests = "$(rpad("Lagrangian Systems",80))" begin include("lagrangian_tests.jl") end
