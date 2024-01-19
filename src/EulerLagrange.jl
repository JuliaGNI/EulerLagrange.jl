module EulerLagrange

    using GeometricBase: NullParameters, StateVariable
    using LinearAlgebra
    using RuntimeGeneratedFunctions
    using Symbolics
    using Symbolics: ArrayLike, FnType, Sym

    import GeometricBase: equations, functions, parameters
    import GeometricEquations: HODE, HODEProblem
    import GeometricEquations: LODE, LODEProblem
    import GeometricEquations: _lode_default_v̄

    RuntimeGeneratedFunctions.init(@__MODULE__)


    export HODE, HODEProblem
    export LODE, LODEProblem

    export equations, functions, parameters, variables


    export symbolize

    include("symbolics.jl")


    export HamiltonianSystem
    export hamiltonian, hamiltonian_variables

    include("hamiltonian.jl")


    export LagrangianSystem, DegenerateLagrangianSystem
    export lagrangian, lagrangian_variables

    include("lagrangian_common.jl")
    include("lagrangian.jl")
    include("lagrangian_degenerate.jl")
    
end
