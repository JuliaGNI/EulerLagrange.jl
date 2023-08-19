using EulerLagrange
using GeometricEquations
using LinearAlgebra
using ModelingToolkit
using Test


L(t,x,v,params) = v ⋅ v / 2 - x ⋅ x / 2

t₀, q₀, v₀ = (0.0, [1.0, 1.0], [0.5, 2.0])
p₀ = zero(v₀)
params = nothing

t, x, v = lagrangian_variables(2)
sym_lag = L(t, x, v, params)
lag_sys = LagrangianSystem(sym_lag, t, x, v)

@test isequal(lagrangian(lag_sys), sym_lag)
@test isequal(variables(lag_sys), (t, x, v))
@test isequal(EulerLagrange.parameters(lag_sys), NamedTuple())

for k in (:L, :EL, :ϑ, :θ, :f, :g, :ω, :Ω, :ϕ, :ψ, :M, :N)
    @test k ∈ keys(EulerLagrange.equations(lag_sys))
end

for k in (:L, :EL, :a, :ϑ, :θ, :f, :g, :ω, :Ω, :ϕ, :ψ, :v̄, :f̄, :M, :P)
    @test k ∈ keys(EulerLagrange.functions(lag_sys))
end


p̃(p, t, q, q̇, params) = p .= q̇
f̃(f, t, q, q̇, params) = f .= -q

p₁, p₂ = zero(p₀), zero(p₀)
f₁, f₂ = zero(p₀), zero(p₀)
    
eqs = functions(lag_sys)

eqs.ϑ(p₁, t₀, q₀, v₀, params)
eqs.f(f₁, t₀, q₀, v₀, params)

p̃(p₂, t₀, q₀, v₀, params)
f̃(f₂, t₀, q₀, v₀, params)

@test eqs.L(t₀, q₀, v₀, params) == L(t₀, q₀, v₀, params)
@test p₁ == p₂
@test f₁ == f₂


ntime = 1000
tstep = 0.01
tspan = (0.0, ntime * tstep)
ics   = (q = q₀, p = p₀, λ = zero(q₀))

lode = LODE(lag_sys)

lprob1 = LODEProblem(lag_sys, tspan, tstep, ics)
lprob2 = LODEProblem(lag_sys, tspan, tstep, q₀, p₀)

@test lode == equation(lprob1) == equation(lprob2)
