```@meta
CurrentModule = EulerLagrange
```

# Euler-Lagrange

This package generates code for the Euler-Lagrange equations as well as Hamilton's equations for [GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl) and related packages.


## Installation

*EulerLagrange.jl* and all of its dependencies can be installed via the Julia REPL by typing 
```
]add EulerLagrange
```

## Usage

Using EulerLagrange.jl is very simple and typically consists of four to five steps:

1) Obtain symbolic variables for a Lagrangian or Hamiltonian system of a given dimension.
2) Obtain a symbolic representation of the parameters of the system if it has any.
3) Build the Lagrangian or Hamiltonian using those symbolic variables and parameters.
4) Construct a `LagrangianSystem` or `HamiltonianSystem`, which is where the actual code generation happens.
5) Generate a `LODEProblem` or `HODEProblem` that can then be solved with [GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl).

Details for the specific system types can be found on the following pages:

```@contents
Pages = [
    "hamiltonian.md",
    "lagrangian.md",
    "degenerate_lagrangian.md",
    "caveats.md",
]
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
