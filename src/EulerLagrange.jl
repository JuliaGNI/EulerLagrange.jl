module EulerLagrange

    using GeometricBase: State
    # using Parameters
    using RuntimeGeneratedFunctions
    using Symbolics
    using Symbolics: FnType, Sym

    import GeometricEquations: HODE, HODEProblem
    import GeometricEquations: LODE, LODEProblem


    include("symbolics.jl")


    export hamiltonian_variables

    include("hamiltonian.jl")


    export LagrangianSystem
    export lagrangian_variables

    include("lagrangian.jl")

    
end
