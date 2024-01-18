```@meta
CurrentModule = EulerLagrange
```

# Lagrangian Systems

The Euler-Lagrange equations, that is the dynamical equations of a Lagrangian system, are given in terms of the Lagrangian ``L(x,v)`` by
```math
\frac{d}{dt} \frac{\partial L}{\partial v} - \frac{\partial L}{\partial x} .
```
For regular (i.e. non-degenerate) Lagrangians, this is a set of second-order ordinary differential equations.
In many numerical applications, it is advantageous to solve the implicit form of these equations, given by
```math
\begin{align}
\frac{d \vartheta}{dt} &= f , &
\vartheta &= \frac{\partial L}{\partial v} , &
f = \frac{\partial L}{\partial x} .
\end{align}
```

In the following, we show how these equations can be obtained for the example of a particle in a square potential.


## Particle in a potential

Before any use, we need to load `EulerLagrange`:
```@example lag
using EulerLagrange
```

Next, we generate symbolic variables for a one-dimensional system:
```@example lag
t, x, v = lagrangian_variables(2)
```

With those variables, we can construct a Lagrangian
```@example lag
using LinearAlgebra
L = v ⋅ v / 2 - x ⋅ x / 2
```

This Lagrangian together with the symbolic variables is then used to construct a `LagrangianSystem`:
```@example lag
lag_sys = LagrangianSystem(L, t, x, v)
```

The constructor computes the Euler-Lagrange equations and generates the corresponding Julia code.
In the last step, we can now construct a `LODEProblem` from the `LagrangianSystem` and some appropriate initial conditions, a time span to integrate over and a time step:
```@example lag
tspan = (0.0, 10.0)
tstep = 0.01

q₀ = [1.0, 1.0]
p₀ = [0.5, 2.0]

lprob = LODEProblem(lag_sys, tspan, tstep, q₀, p₀)
```

Should we fancy so, we can integrate this system using GeometricIntegrators:
```@example lag
using GeometricIntegrators
sol = integrate(lprob, Gauss(1))

using CairoMakie
fig = lines(parent(sol.q[:,1]), parent(sol.q[:,2]);
    axis = (; xlabel = "x₁", ylabel = "x₂", title = "Particle moving in a square potential"),
    figure = (; size = (800,600), fontsize = 22))
save("particle_vi.svg", fig); nothing # hide
```

![](particle_vi.svg)


```@setup lag_params
using CairoMakie
using EulerLagrange
using LinearAlgebra
using GeometricIntegrators

t, x, v = lagrangian_variables(2)

tspan = (0.0, 10.0)
tstep = 0.01

q₀ = [1.0, 1.0]
p₀ = [0.5, 2.0]
```

We can also include parametric dependencies in the Lagrangian.
Consider, for example, a parameter `α` that determines the strength of the potential.

The easiest way, to account for parameters, is to create a named tuple with typical values for each parameters, e.g.,
```@example lag_params
params = (α = 5.0,)
```

In the next step, we use the function `symbolize` to generate a symbolic version of the parameters:
```@example lag_params
sparams = symbolize(params)
```

Now we modify the Lagrangian to account for the parameter:
```@example lag_params
L = v ⋅ v / 2 - sparams.α * (x ⋅ x) / 2
```

From here on, everything follows along the same lines as before, the only difference being that we also need to pass the symbolic parameters `sparams` to the `LagrangianSystem` constructor:
```@example lag_params
lag_sys = EulerLagrange.LagrangianSystem(L, t, x, v, sparams)
```

Analogously, we need to pass actual parameter values `params` to the `LODEProblem` constructor via the `parameters` keyword argument:
```@example lag_params
lprob = LODEProblem(lag_sys, tspan, tstep, q₀, p₀; parameters = params)
```

This problem can again be integrated using GeometricIntegrators:
```@example lag_params
sol = integrate(lprob, Gauss(1))

fig = lines(parent(sol.q[:,1]), parent(sol.q[:,2]);
    axis = (; xlabel = "x₁", ylabel = "x₂", title = "Particle moving in a square potential"),
    figure = (; size = (800,600), fontsize = 22))
save("particle_vi_param.svg", fig); nothing # hide
```

![](particle_vi_param.svg)
