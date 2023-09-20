using EulerLagrange
using GeometricEquations
using LinearAlgebra
using Test


H(t, x, v, params) = params.a₁ * x[1] + params.a₂ * x[2] + params.b₁ * log(x[1]) + params.b₂ * log(x[2])
L(t, x, v, params) = log(x[2]) / x[1] / 2 * v[1] - log(x[1]) / x[2] / 2 * v[2] - H(t,x,v,params)

t₀, q₀ = 0.0, [2.0, 1.0]

params = (
    a₁ = -1.0,
    a₂ = -1.0,
    b₁ = 1.0,
    b₂ = 2.0,
)

t, x, v = lagrangian_variables(2)
sparams = symbolize(params)
lag_sys = LagrangianSystem(L(t,x,v,sparams), t, x, v, sparams)

dHd₁(t, q, params) = params.a₁ + params.b₁ / q[1]
dHd₂(t, q, params) = params.a₂ + params.b₂ / q[2]
dϑ₁dx₁(t, q) = - log(q[2]) / q[1]^2 / 2
dϑ₁dx₂(t, q) = + 1 / (q[1] * q[2]) / 2
dϑ₂dx₁(t, q) = - 1 / (q[2] * q[1]) / 2
dϑ₂dx₂(t, q) = + log(q[1]) / q[2]^2 / 2
f₁(t, q, v) = dϑ₁dx₁(t,q) * v[1] + dϑ₂dx₁(t,q) * v[2]
f₂(t, q, v) = dϑ₁dx₂(t,q) * v[1] + dϑ₂dx₂(t,q) * v[2]

function p̃(p, t, q, q̇, params)
    p[1] = + log(q[2]) / q[1] / 2
    p[2] = - log(q[1]) / q[2] / 2
end

function f̃(f, t, q, q̇, params)
    f[1] = f₁(t,q,q̇) - dHd₁(t, q, params)
    f[2] = f₂(t,q,q̇) - dHd₂(t, q, params)
end

t₀, q₀, v₀ = (0.0, [1.0, 1.0], [0.5, 2.0])
p₀ = zero(v₀)

p₁, p₂ = zero(p₀), zero(p₀)
ṗ₁, ṗ₂ = zero(p₀), zero(p₀)
    
eqs = functions(lag_sys)

eqs.ϑ(p₁, t₀, q₀, v₀, params)
eqs.f(ṗ₁, t₀, q₀, v₀, params)

p̃(p₂, t₀, q₀, v₀, params)
f̃(ṗ₂, t₀, q₀, v₀, params)

@test eqs.L(t₀, q₀, v₀, params) == L(t₀, q₀, v₀, params)
@test p₁ == p₂
@test ṗ₁ == ṗ₂
