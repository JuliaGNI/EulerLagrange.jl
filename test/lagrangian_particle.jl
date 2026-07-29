using EulerLagrange
using GeometricEquations
using LinearAlgebra
using Symbolics
using Test

# Values produced by the generated functions are compared against independently written reference
# expressions with a relative tolerance, not with `==`. The two evaluate the same mathematical
# quantity in a different order, and `cse=true` additionally binds constant nodes to temporaries,
# so results may differ in the last bit.
const RTOL = 1.0e-14


# Initial conditions and general variables

t, x, v = lagrangian_variables(2)

t₀, q₀, v₀ = (0.0, [1.0, 1.0], [0.5, 2.0])

p₀ = zero(v₀)
p₁, p₂ = zero(p₀), zero(p₀)
f₁, f₂ = zero(p₀), zero(p₀)

ntime = 1000
tstep = 0.01
tspan = (0.0, ntime * tstep)
ics = (q=StateVariable(q₀), p=StateVariable(p₀), λ=AlgebraicVariable(zero(q₀)))


# Test without parameters

L(t, x, v, params) = v ⋅ v / 2 - x ⋅ x / 2

params = nothing

sym_lag = L(t, x, v, params)
lag_sys = LagrangianSystem(sym_lag, t, x, v)

# @test isequal(lagrangian(lag_sys), simplify(sym_lag)) # TODO: reactivate!
@test isequal(variables(lag_sys), (t, x, v))
@test isequal(EulerLagrange.parameters(lag_sys), NamedTuple())

for k in (:L, :EL, :f, :g, :ϑ, :ω, :ϕ, :ψ, :M, :N)
    @test k ∈ keys(EulerLagrange.equations(lag_sys))
end

for k in (:L, :EL, :f, :g, :p, :ϑ, :ω, :ϕ, :ψ, :M)# :a, :P
    @test k ∈ keys(EulerLagrange.functions(lag_sys))
end


p̃(p, t, q, q̇, params) = p .= q̇
f̃(f, t, q, q̇, params) = f .= -q

eqs = functions(lag_sys)

eqs.ϑ(p₁, t₀, q₀, v₀, params)
eqs.f(f₁, t₀, q₀, v₀, params)

p̃(p₂, t₀, q₀, v₀, params)
f̃(f₂, t₀, q₀, v₀, params)

@test eqs.L(t₀, q₀, v₀, params) ≈ L(t₀, q₀, v₀, params) rtol = RTOL
@test p₁ ≈ p₂ rtol = RTOL
@test f₁ ≈ f₂ rtol = RTOL


lode = LODE(lag_sys)

lprob1 = LODEProblem(lag_sys, tspan, tstep, ics)
lprob2 = LODEProblem(lag_sys, tspan, tstep, q₀, p₀)

@test lode == equation(lprob1) == equation(lprob2)


# Test with parameters

Lₚ(t, x, v, params) = v ⋅ v / 2 - params.α * (x ⋅ x) / 2

params = (α=5.0,)
sparams = symbolize(params)

sym_lag = Lₚ(t, x, v, params)
lag_sys = LagrangianSystem(sym_lag, t, x, v, sparams)

# @test isequal(lagrangian(lag_sys), sym_lag) # TODO: reactivate!
@test isequal(EulerLagrange.parameters(lag_sys), sparams)

p̃ₚ(p, t, q, q̇, params) = p .= q̇
f̃ₚ(f, t, q, q̇, params) = f .= -params.α * q

eqs = functions(lag_sys)

eqs.ϑ(p₁, t₀, q₀, v₀, params)
eqs.f(f₁, t₀, q₀, v₀, params)

p̃ₚ(p₂, t₀, q₀, v₀, params)
f̃ₚ(f₂, t₀, q₀, v₀, params)

@test eqs.L(t₀, q₀, v₀, params) ≈ Lₚ(t₀, q₀, v₀, params) rtol = RTOL
@test p₁ ≈ p₂ rtol = RTOL
@test f₁ ≈ f₂ rtol = RTOL
