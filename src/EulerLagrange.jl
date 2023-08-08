module EulerLagrange

    # using Parameters
    using RuntimeGeneratedFunctions
    using Symbolics
    using Symbolics: Sym


    include("symbolics.jl")

    export HamiltonianSystem
    export hamiltonian_variables

    include("hamiltonian.jl")


    export LagrangianSystem
    export lagrangian_variables

    include("lagrangian.jl")

    
end
