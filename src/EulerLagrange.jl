module EulerLagrange

using GeometricBase: NullParameters, StateVariable, AlgebraicVariable
using LinearAlgebra
using RuntimeGeneratedFunctions
using Symbolics
using Symbolics: FnType, Sym

import GeometricBase: equations, functions, parameters
import GeometricEquations: HODE, HODEProblem
import GeometricEquations: LODE, LODEProblem, LDAE, LDAEProblem
import GeometricEquations: ODE, ODEProblem
import GeometricEquations: _lode_default_v̄

RuntimeGeneratedFunctions.init(@__MODULE__)


export ODE, ODEProblem
export HODE, HODEProblem
export LODE, LODEProblem
export LDAE, LDAEProblem

export equations, functions, parameters, variables


export symbolize

include("symbolics.jl")


export HamiltonianSystem
export hamiltonian, hamiltonian_variables, hamiltonian_derivatives

include("hamiltonian.jl")


export LagrangianSystem, DegenerateLagrangianSystem
export lagrangian, lagrangian_variables, lagrangian_derivatives

include("lagrangian_common.jl")
include("lagrangian.jl")
include("lagrangian_degenerate.jl")

end
