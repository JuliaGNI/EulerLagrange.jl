using EulerLagrange
using GeometricEquations
using LinearAlgebra
using Symbolics
using Test


H(t, q, p, params) = p ⋅ p / 2 + q ⋅ q / 2

t₀, q₀, p₀ = 0.0, [1.0, 1.0], [0.5, 2.0]
z₀ = [q₀..., p₀...]
params = nothing

t, q, p = hamiltonian_variables(2)
sym_ham = H(t, q, p, params)
ham_sys = HamiltonianSystem(sym_ham, t, q, p)

@test isequal(hamiltonian(ham_sys), Symbolics.simplify(sym_ham))
@test isequal(variables(ham_sys), (t, q, p))
@test isequal(EulerLagrange.parameters(ham_sys), NamedTuple())

for k in (:H, :EH, :EHq, :EHp, :v, :f, :ż)
    @test k ∈ keys(EulerLagrange.equations(ham_sys))
    @test k ∈ keys(EulerLagrange.functions(ham_sys))
end

ṽ(v, t, q, p, params) = v .= p
f̃(f, t, q, p, params) = f .= -q
ż̃(ż, t, q, p, params) = ż .= [p..., -q...]

v₁, v₂ = zero(q₀), zero(q₀)
f₁, f₂ = zero(p₀), zero(p₀)
ż₁, ż₂ = zero(z₀), zero(z₀)
    
eqs = functions(ham_sys)

eqs.v(v₁, t₀, q₀, p₀, params)
eqs.f(f₁, t₀, q₀, p₀, params)
eqs.ż(ż₁, t₀, q₀, p₀, params)

ṽ(v₂, t₀, q₀, p₀, params)
f̃(f₂, t₀, q₀, p₀, params)
ż̃(ż₂, t₀, q₀, p₀, params)
    
@test eqs.H(t₀, q₀, p₀, params) == H(t₀, q₀, p₀, params)
@test v₁ == v₂
@test f₁ == f₂
@test ż₁ == ż₂


ntime = 1000
tstep = 0.01
tspan = (0.0, ntime * tstep)
ics   = (q = StateVariable(q₀), p = StateVariable(p₀))

hode = HODE(ham_sys)

hprob1 = HODEProblem(ham_sys, tspan, tstep, ics)
hprob2 = HODEProblem(ham_sys, tspan, tstep, q₀, p₀)

@test hode == equation(hprob1) == equation(hprob2)
