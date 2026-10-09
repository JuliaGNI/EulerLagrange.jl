using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Symbolize" include("symbolics.jl")
    @safetestset "Hamiltonian: general functionality" include("hamiltonian.jl")
    @safetestset "Hamiltonian: particle in square potential" include("integration/hamiltonian_particle.jl")
    @safetestset "Hamiltonian: harmonic oscillator" include("integration/hamiltonian_oscillator.jl")
    @safetestset "Lagrangian: general functionality" include("lagrangian_common.jl")
    @safetestset "Lagrangian: particle in square potential" include("integration/lagrangian_particle.jl")
    @safetestset "Lagrangian: Lotka-Volterra" include("integration/lagrangian_lotka_volterra.jl")
    @safetestset "Lagrangian: solar system" include("integration/lagrangian_solar_system.jl")
    @safetestset "Common subexpression elimination" include("integration/cse_tests.jl")
    @safetestset "Generated function signatures" include("integration/signature_tests.jl")
end
