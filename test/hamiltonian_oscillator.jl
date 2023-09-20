using EulerLagrange
using GeometricEquations
using LinearAlgebra
using Test


H(t, q, p, params) = dot(p,p) / 2 + params.k * dot(q,q) / 2

t₀, q₀, p₀ = 0.0, [0.5], [0.0]
params = (k=0.5, ω=√0.5)

t, q, p = hamiltonian_variables(1)
sparams = symbolize(params)
ham_sys = HamiltonianSystem(H(t, q, p, sparams), t, q, p, sparams)

ṽ(v, t, q, p, params) = v .= p
f̃(f, t, q, p, params) = f .= - params.k .* q

v₁, v₂ = zero(q₀), zero(q₀)
f₁, f₂ = zero(p₀), zero(p₀)
    
eqs = functions(ham_sys)

eqs.v(v₁, t₀, q₀, p₀, params)
eqs.f(f₁, t₀, q₀, p₀, params)

ṽ(v₂, t₀, q₀, p₀, params)
f̃(f₂, t₀, q₀, p₀, params)
    
@test eqs.H(t₀, q₀, p₀, params) == H(t₀, q₀, p₀, params)
@test v₁ == v₂
@test f₁ == f₂
