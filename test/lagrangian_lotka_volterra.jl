using EulerLagrange
using LinearAlgebra
using Test

# Values produced by the generated functions are compared against independently written reference
# expressions with a relative tolerance, not with `==`. The two evaluate the same mathematical
# quantity in a different order, and `cse=true` additionally binds constant nodes to temporaries,
# so results may differ in the last bit.
const RTOL = 1.0e-14


ϑ(t, x, v, params) = [log(x[2]) / x[1] / 2, -log(x[1]) / x[2] / 2]
H(t, x, v, params) = params.a₁ * x[1] + params.a₂ * x[2] + params.b₁ * log(x[1]) + params.b₂ * log(x[2])
K(t, x, v, params) = ϑ(t, x, v, params) ⋅ v
L(t, x, v, params) = K(t, x, v, params) - H(t, x, v, params)

dHd₁(t, q, params) = params.a₁ + params.b₁ / q[1]
dHd₂(t, q, params) = params.a₂ + params.b₂ / q[2]
dϑ₁dx₁(t, q) = -log(q[2]) / q[1]^2 / 2
dϑ₁dx₂(t, q) = +1 / (q[1] * q[2]) / 2
dϑ₂dx₁(t, q) = -1 / (q[2] * q[1]) / 2
dϑ₂dx₂(t, q) = +log(q[1]) / q[2]^2 / 2


v₁(t, q, params) = +q[1] * (params.a₂ * q[2] + params.b₂)
v₂(t, q, params) = -q[2] * (params.a₁ * q[1] + params.b₁)
f₁(t, q, v) = dϑ₁dx₁(t, q) * v[1] + dϑ₂dx₁(t, q) * v[2]
f₂(t, q, v) = dϑ₁dx₂(t, q) * v[1] + dϑ₂dx₂(t, q) * v[2]

function ṽ(v, t, q, params)
    v[1] = v₁(t, q, params)
    v[2] = v₂(t, q, params)
end

function p̃(p, t, q, q̇, params)
    p[1] = +log(q[2]) / q[1] / 2
    p[2] = -log(q[1]) / q[2] / 2
end

function f̃(f, t, q, q̇, params)
    f[1] = f₁(t, q, q̇) - dHd₁(t, q, params)
    f[2] = f₂(t, q, q̇) - dHd₂(t, q, params)
end

t₀, q₀, v₀ = (0.0, [2.0, 1.0], [0.5, 2.0])
p₀ = zero(v₀)
λ₀ = zero(v₀)

tspan = (0.0, 1.0)
tstep = 0.1

params = (
    a₁=-1.0,
    a₂=-1.0,
    b₁=1.0,
    b₂=2.0,
)

params_alt = (
    a₁=-1.0,
    a₂=-1.0,
    b₁=2.0,
    b₂=1.0,
)

p̃(p₀, t₀, q₀, v₀, params)


# Symbolic variables and parameters

t, x, v = lagrangian_variables(2)
sparams = symbolize(params)


# LagrangianSystem

lag_sys = LagrangianSystem(L(t, x, v, sparams), t, x, v, sparams)

p₁, p₂ = zero(p₀), zero(p₀)
ṗ₁, ṗ₂ = zero(p₀), zero(p₀)

eqs = functions(lag_sys)

eqs.ϑ(p₁, t₀, q₀, v₀, params)
eqs.f(ṗ₁, t₀, q₀, v₀, params)

p̃(p₂, t₀, q₀, v₀, params)
f̃(ṗ₂, t₀, q₀, v₀, params)

@test eqs.L(t₀, q₀, v₀, params) ≈ L(t₀, q₀, v₀, params) rtol = RTOL
@test p₁ ≈ p₂ atol = 2eps()
@test ṗ₁ ≈ ṗ₂ atol = 2eps()


@test eqs.L(t₀, q₀, v₀, params_alt) != L(t₀, q₀, v₀, params)

@test_nowarn LODE(lag_sys)
@test_nowarn LODEProblem(lag_sys, tspan, tstep, q₀, v₀; parameters=params)


# DegenerateLagrangianSystem

deg_lag_sys = DegenerateLagrangianSystem(K(t, x, v, sparams), H(t, x, v, sparams), t, x, v, sparams)

q̇₁, q̇₂ = zero(q₀), zero(q₀)
p₁, p₂ = zero(p₀), zero(p₀)
ṗ₁, ṗ₂ = zero(p₀), zero(p₀)

deg_eqs = functions(deg_lag_sys)

deg_eqs.ẋ(q̇₁, t₀, q₀, params)
deg_eqs.ϑ(p₁, t₀, q₀, v₀, params)
deg_eqs.f(ṗ₁, t₀, q₀, v₀, params)

ṽ(q̇₂, t₀, q₀, params)
p̃(p₂, t₀, q₀, v₀, params)
f̃(ṗ₂, t₀, q₀, v₀, params)

@test deg_eqs.L(t₀, q₀, v₀, params) ≈ L(t₀, q₀, v₀, params) rtol = RTOL
@test q̇₁ ≈ q̇₂ atol = 2eps()
@test p₁ ≈ p₂ atol = 2eps()
@test ṗ₁ ≈ ṗ₂ atol = 2eps()

@test deg_eqs.L(t₀, q₀, v₀, params_alt) != L(t₀, q₀, v₀, params)


# The secondary constraint ψ = ṗ - q̇ ⋅ ∇ϑ, called the way `LDAE` calls it:
# ψ(out, t, q, v, p, q̇, ṗ, params). The algebraic velocity `v` and the time derivative `q̇` are
# separate arguments, so they are chosen distinct here — contracting ∇ϑ with the wrong one is
# exactly the failure this pins down.

function ψ̃(ψ, t, q, q̇, ṗ)
    ψ[1] = ṗ[1] - (dϑ₁dx₁(t, q) * q̇[1] + dϑ₁dx₂(t, q) * q̇[2])
    ψ[2] = ṗ[2] - (dϑ₂dx₁(t, q) * q̇[1] + dϑ₂dx₂(t, q) * q̇[2])
end

q̇₀ = [0.7, -0.3]
ṗ₀ = [1.3, 0.9]

ψ₁, ψ₂ = zero(p₀), zero(p₀)

deg_eqs.ψ(ψ₁, t₀, q₀, v₀, p₀, q̇₀, ṗ₀, params)
ψ̃(ψ₂, t₀, q₀, q̇₀, ṗ₀)

@test ψ₁ ≈ ψ₂ rtol = RTOL

# ψ must depend on q̇, not on the algebraic velocity v
ψ₃ = zero(p₀)
deg_eqs.ψ(ψ₃, t₀, q₀, 2 .* v₀, p₀, q̇₀, ṗ₀, params)
@test ψ₃ == ψ₁


@test_nowarn ODE(deg_lag_sys)
@test_nowarn ODEProblem(deg_lag_sys, tspan, tstep, q₀; parameters=params)

@test_nowarn LODE(deg_lag_sys)
@test_nowarn LODEProblem(deg_lag_sys, tspan, tstep, q₀, p₀; parameters=params)

@test_nowarn LDAE(deg_lag_sys)
@test_nowarn LDAEProblem(deg_lag_sys, tspan, tstep, q₀, p₀, λ₀; parameters=params)
