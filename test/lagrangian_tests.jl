@safetestset LagrangianGeneral = "$(rpad("General Functionality",80))" begin include("lagrangian_general.jl") end
@safetestset LagrangianParticle = "$(rpad("Particle in square potential",80))" begin include("lagrangian_particle.jl") end
@safetestset LagrangianLotkaVolterra = "$(rpad("Lotka-Volterra",80))" begin include("lagrangian_lotka_volterra.jl") end
