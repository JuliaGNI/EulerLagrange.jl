@safetestset HamiltonianGeneral = "$(rpad("General Functionality",80))" begin include("hamiltonian_general.jl") end
@safetestset HamiltonianParticle = "$(rpad("Particle in square potential",80))" begin include("hamiltonian_particle.jl") end
@safetestset HamiltonianOscillator = "$(rpad("Harmonic oscillator",80))" begin include("hamiltonian_oscillator.jl") end
