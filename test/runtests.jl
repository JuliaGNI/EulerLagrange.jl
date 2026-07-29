using EulerLagrange
using SafeTestsets

@safetestset SymbolizeTests = "$(rpad("Symbolize",80))" begin
    include("symbolize_tests.jl")
end
@safetestset HamiltonianTests = "$(rpad("Hamiltonian Systems",80))" begin
    include("hamiltonian_tests.jl")
end
@safetestset LagrangianTests = "$(rpad("Lagrangian Systems",80))" begin
    include("lagrangian_tests.jl")
end
@safetestset CSETests = "$(rpad("Common Subexpression Elimination",80))" begin
    include("cse_tests.jl")
end
