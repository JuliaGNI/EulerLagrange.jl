module EulerLagrange

    using GeometricBase: NullParameters, State
    using RuntimeGeneratedFunctions
    using Symbolics
    using Symbolics: FnType, Sym

    import GeometricBase: equations, functions, parameters
    import GeometricEquations: HODE, HODEProblem
    import GeometricEquations: LODE, LODEProblem


    export equations, functions, parameters, variables


    export symbolize

    include("symbolics.jl")


    export HamiltonianSystem
    export hamiltonian, hamiltonian_variables

    include("hamiltonian.jl")


    export LagrangianSystem
    export lagrangian, lagrangian_variables

    include("lagrangian.jl")
    
end
