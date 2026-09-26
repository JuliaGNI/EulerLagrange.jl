using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Symbolize" include("symbolics.jl")
    @safetestset "Hamiltonian: general functionality" include("hamiltonian_general.jl")
    @safetestset "Hamiltonian: particle in square potential" include("hamiltonian_particle.jl")
    @safetestset "Hamiltonian: harmonic oscillator" include("hamiltonian_oscillator.jl")
    @safetestset "Lagrangian: general functionality" include("lagrangian_general.jl")
    @safetestset "Lagrangian: particle in square potential" include("lagrangian_particle.jl")
    @safetestset "Lagrangian: Lotka-Volterra" include("lagrangian_lotka_volterra.jl")
    @safetestset "Lagrangian: solar system" include("lagrangian_solar_system.jl")
    @safetestset "Common subexpression elimination" include("cse_tests.jl")
    @safetestset "Generated function signatures" include("signature_tests.jl")
end
