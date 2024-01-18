# Hamiltonian Systems

Hamilton's equations of motion are given in terms of the Hamiltonian ``H(q,p)`` by
```math
\begin{align*}
\frac{dq}{dt} &= \frac{\partial H}{\partial p} , &
\frac{dp}{dt} &= - \frac{\partial H}{\partial q} .
\end{align*}
```

In the following, we show how these equations can be obtained for the example of a harmonic oscillator.


## Harmonic Oscillator

Before any use, we need to load `EulerLagrange`:
```@example ham
using EulerLagrange
```

Next, we generate symbolic variables for a one-dimensional system:
```@example ham
t, q, p = hamiltonian_variables(1)
```

We define a named tuple with typical values for the parameters, e.g.,
```@example ham
params = (k=0.5, ω=√0.5)
```

We use the function `symbolize` to generate a symbolic version of the parameters:
```@example ham
sparams = symbolize(params)
```

Now we can define the Hamiltonian function:
```@example ham
using LinearAlgebra
H(t, q, p, params) = p ⋅ p / 2 + params.k * (q ⋅ q) / 2
```

The Hamiltonian, evaluated on and together with the symbolic variables and parameters is used to construct a `HamiltonianSystem`:
```@example ham
ham_sys = HamiltonianSystem(H(t, q, p, sparams), t, q, p, sparams)
```
The constructor computes Hamilton's equations and generates the corresponding Julia code.
In the last step, we can now construct a `HODEProblem` from the `HamiltonianSystem` and some appropriate initial conditions, a time span to integrate over and a time step:
```@example ham
tspan = (0.0, 10.0)
tstep = 0.01
q₀, p₀ = [0.5], [0.0]

hprob = HODEProblem(ham_sys, tspan, tstep, q₀, p₀; parameters = params)
```

We can integrate this system using GeometricIntegrators:
```@example ham
using GeometricIntegrators
sol = integrate(hprob, Gauss(1))

using CairoMakie
fig = lines(parent(sol.q[:,1]), parent(sol.p[:,1]);
    axis = (; xlabel = "q₁", ylabel = "p₁", title = "Harmonic Oscillator"),
    figure = (; size = (800,600), fontsize = 22))
save("harmonic-oscillator-hamiltonian.svg", fig); nothing # hide
```

![](harmonic-oscillator-hamiltonian.svg)
