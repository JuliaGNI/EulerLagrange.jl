module EulerLagrange

    using GeometricBase: State
    # using Parameters
    using RuntimeGeneratedFunctions
    using Symbolics
    using Symbolics: FnType, Sym

    import GeometricBase: equation, equations
    import GeometricEquations: HODE, HODEProblem
    import GeometricEquations: LODE, LODEProblem


    export equation, equations


    include("symbolics.jl")


    export HamiltonianSystem
    export hamiltonian_variables

    include("hamiltonian.jl")


    export LagrangianSystem
    export lagrangian_variables

    include("lagrangian.jl")
    
end
