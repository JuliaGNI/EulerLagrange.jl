# EulerLagrange

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNI.github.io/EulerLagrange.jl/stable)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaGNI.github.io/EulerLagrange.jl/latest)
[![Build Status](https://github.com/JuliaGNI/EulerLagrange.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/EulerLagrange.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGNI/EulerLagrange.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/EulerLagrange.jl)

This package generates code for the Euler-Lagrange equations as well as Hamilton's equations for [GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl) and related packages.


## Installation

*EulerLagrange.jl* and all of its dependencies can be installed via the Julia REPL by typing 
```
]add EulerLagrange
```

## Basic usage

Using EulerLagrange.jl is very simple and typically consists of four to five steps:

1) Obtain symbolic variables for a Lagrangian or Hamiltonian system of a given dimension.
2) Obtain a symbolic representation of the parameters of the system if it has any.
3) Build the Lagrangian or Hamiltonian using those symbolic variables and parameters.
4) Construct a `LagrangianSystem` or `HamiltonianSystem`, which is where the actual code generation happens.
5) Generate a `LODEProblem` or `HODEProblem` that can then be solved with [GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl).

In the following, we showcase this procedure for a particle in a square potential.

Before any use, we need to load `EulerLagrange`:
```julia
using EulerLagrange
```

Next, we generate symbolic variables for a one-dimensional system:
```julia
t, x, v = lagrangian_variables(1)
```

With those variables, we can construct a Lagrangian
```julia
L = v ⋅ v / 2 - x ⋅ x / 2
```

This Lagrangian together with the symbolic variables is then used to construct a `LagrangianSystem`:
```julia
lag_sys = LagrangianSystem(L, t, x, v)
```

The constructor computes the Euler-Lagrange equations and generates the corresponding Julia code.
In the last step, we can now construct a `LODEProblem` from the `LagrangianSystem` and some appropriate initial conditions, a time span to integrate over and a time step:
```julia
tspan = (0.0, 10.0)
tstep = 0.01

q₀ = [1.0]
p₀ = [0.5]

lprob = LODEProblem(lag_sys, tspan, tstep, q₀, p₀)
```

Should we fancy so, we can integrate this system using GeometricIntegrators:
```julia
using GeometricIntegrators
integrate(lprob, ExplicitMidpoint())
```


## References

If you use EulerLagrange.jl in your work, please consider citing it by

```
@misc{Kraus:2023:EulerLagrange,
  title={EulerLagrange.jl: Code generation for Euler-Lagrange equations in Julia},
  author={Kraus, Michael},
  year={2023},
  howpublished={\url{https://github.com/JuliaGNI/EulerLagrange.jl}},
  doi={10.5281/zenodo.8241048}
}
```
