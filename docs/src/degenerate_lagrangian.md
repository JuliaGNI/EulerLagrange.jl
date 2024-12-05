# Degenerate Lagrangian Systems

A Lagrangian ``L(x,v)`` is said to be degenerate if the determinant of the Hessian with respect to the velocities vanishes, i.e.,
```math
\left\vert \frac{\partial^2 L}{\partial v^i \, \partial v^j} \right\vert = 0 .
```
The Euler-Lagrange equations,
```math
\frac{d}{dt} \frac{\partial L}{\partial v} - \frac{\partial L}{\partial x} = 0 ,
```
of such Lagrangians are a set of first-order ordinary differential equations and not second-order differential equations as for non-degenerate Lagrangians.

Of particular interest in various applications are degenerate Lagrangians of the form
```math
L(x,v) = \vartheta(x) \cdot v - H (x) .
```
Their Euler-Lagrange equations take the form
```math
\frac{d \vartheta}{dt} (x) = \nabla \vartheta(x) \cdot v - \nabla H (x) .
```

We exemplify this with the Lotka-Volterra problem in 2d.


## Lotka-Volterra

Before any use, we need to load `EulerLagrange`:
```@example deglag
using EulerLagrange
using LinearAlgebra
```

Next, we generate symbolic variables for a two-dimensional Lagrangian system:
```@example deglag
t, x, v = lagrangian_variables(2)
```

We define a named tuple with typical values for the parameters, e.g.,
```@example deglag
params = (
    a₁ = -1.0,
    a₂ = -1.0,
    b₁ = 1.0,
    b₂ = 2.0,
)
```

We use the function `symbolize` to generate a symbolic version of the parameters:
```@example deglag
sparams = symbolize(params)
```

Define the Hamiltonian function and the symplectic potential:
```@example deglag
ϑ(x, params) = [log(x[2]) / x[1] / 2, - log(x[1]) / x[2] / 2]
H(x, params) = params.a₁ * x[1] + params.a₂ * x[2] + params.b₁ * log(x[1]) + params.b₂ * log(x[2])
```

The Hamiltonian and the symplectic potential, evaluated on and together with the symbolic variables and parameters are used to construct a `DegenerateLagrangianSystem`:
```@example deglag
lag_sys = DegenerateLagrangianSystem(ϑ(x,sparams) ⋅ v, H(x,sparams), t, x, v, sparams)
```
The constructor computes the Euler-Lagrange equations and generates the corresponding Julia code.

In the last step, we can now construct a `LODEProblem` from the `LagrangianSystem` and some appropriate initial conditions, a time span to integrate over and a time step:
```@example deglag
tspan = (0.0, 10.0)
tstep = 0.01

q₀ = [2.0, 1.0]
p₀ = ϑ(q₀, params)

lprob = LODEProblem(lag_sys, tspan, tstep, q₀, p₀; parameters = params)
```

We can integrate this system using GeometricIntegrators:
```@example deglag
using GeometricIntegrators
sol = integrate(lprob, Gauss(1))

using CairoMakie
fig = lines(parent(sol.q[:,1]), parent(sol.q[:,2]);
    axis = (; xlabel = "x₁", ylabel = "x₂", title = "Lotka-Volterra system in 2d"),
    figure = (; size = (800,600), fontsize = 22))
save("lotka-volterra.svg", fig); nothing # hide
```

![](lotka-volterra.svg)
