using EulerLagrange
using LinearAlgebra
using Test

# `cse=true` only names repeated subexpressions; it never reassociates or reorders operands.
# Results are therefore equal up to floating-point rounding — `SymbolicUtils.Code.cse` also binds
# constant nodes to temporaries, which can shift the association by a rounding step — so these
# comparisons use a tolerance rather than `==`.

const RTOL = 1.0e-12

# A gravitational N-body model: every pairwise distance appears in many components of the
# derivatives, so this is a case where CSE has real work to do.
const G = 2.95912208286e-4
const masses = [1.0, 9.54786104043e-4, 2.85583733151e-4, 4.37273164546e-5]

const N = length(masses)
const DIM = 3
const D = N * DIM

function potential(q, params)
    Q = reshape(q, DIM, N)
    s = zero(eltype(q))
    for i in 2:N, j in 1:i-1
        Δ₁ = Q[1, i] - Q[1, j]
        Δ₂ = Q[2, i] - Q[2, j]
        Δ₃ = Q[3, i] - Q[3, j]
        s += params.G * params.m[i] * params.m[j] / sqrt(Δ₁^2 + Δ₂^2 + Δ₃^2)
    end
    -s
end

function kinetic(v, params)
    V = reshape(v, DIM, N)
    s = zero(eltype(v))
    for i in 1:N
        s += params.m[i] * (V[1, i]^2 + V[2, i]^2 + V[3, i]^2)
    end
    s / 2
end

lagrangian(t, q, v, params) = kinetic(v, params) - potential(q, params)

function hamiltonian(t, q, p, params)
    P = reshape(p, DIM, N)
    s = zero(eltype(p))
    for i in 1:N
        s += (P[1, i]^2 + P[2, i]^2 + P[3, i]^2) / params.m[i]
    end
    s / 2 + potential(q, params)
end

const params = (G=G, m=masses)

# Deterministic, non-degenerate states: no two bodies coincide, so the potential is finite.
const t₀ = 0.7
const q₀ = collect(range(-8.0, 11.0; length=D)) .+ 0.25 .* collect(1.0:D)
const v₀ = collect(range(-0.004, 0.006; length=D))
const p₀ = 0.5 .* v₀
const p₁ = collect(range(-0.3, 0.4; length=D))
const a₀ = collect(range(-2.0e-8, 3.0e-8; length=D))  # accelerations, the `Λ` slot

"Evaluate `key` from both variants and return the two results."
function evaluate(withcse, without, key, sizes, args...)
    a = zeros(sizes)
    b = zeros(sizes)
    getproperty(withcse, key)(a, t₀, args..., params)
    getproperty(without, key)(b, t₀, args..., params)
    return a, b
end

"Compare `key` built with `cse=true` against `cse=false`, given the arguments after `t`."
function compare(withcse, without, key, sizes, args...)
    a, b = evaluate(withcse, without, key, sizes, args...)
    @test a ≈ b rtol = RTOL
    # a quantity that evaluates to all zeros would make the comparison above vacuous
    @test any(!iszero, b)
end

@testset "$(rpad("Lagrangian system: cse=true matches cse=false", 80))" begin
    t, x, v = lagrangian_variables(D)
    sparams = symbolize(params)
    sym_lag = lagrangian(t, x, v, sparams)

    with = functions(LagrangianSystem(sym_lag, t, x, v, sparams; simplify=false, cse=true))
    without = functions(LagrangianSystem(sym_lag, t, x, v, sparams; simplify=false, cse=false))

    @test with.L(t₀, q₀, v₀, params) ≈ without.L(t₀, q₀, v₀, params) rtol = RTOL
    @test with.L(t₀, q₀, v₀, params) ≈ lagrangian(t₀, q₀, v₀, params) rtol = RTOL
    @test with.p(t₀, q₀, v₀, params) ≈ without.p(t₀, q₀, v₀, params) rtol = RTOL

    for key in (:f, :ϑ)
        compare(with, without, key, D, q₀, v₀)
    end
    compare(with, without, :ϕ, D, q₀, v₀, p₀)
    compare(with, without, :M, (D, D), q₀, v₀)

    # `g`, `EL = f - g` and `ψ = F - g` take the acceleration in their trailing `Λ` slot.
    compare(with, without, :g, D, q₀, v₀, a₀)
    compare(with, without, :EL, D, q₀, v₀, a₀)
    compare(with, without, :ψ, D, q₀, v₀, p₀, p₁, a₀)

    # `ω = dθ_L` for the Lagrange one-form θ_L = ϑᵢ dqⁱ. A `LagrangianSystem` is a regular
    # (second-order) system, so `ω` is the full 2n×2n two-form on (q, q̇); here it reduces to the
    # canonical `[0 -M; M 0]` with M = diag(m).
    compare(with, without, :ω, (2D, 2D), q₀, v₀)
end

@testset "$(rpad("Hamiltonian system: cse=true matches cse=false", 80))" begin
    t, q, p = hamiltonian_variables(D)
    sparams = symbolize(params)
    sym_ham = hamiltonian(t, q, p, sparams)

    with = functions(HamiltonianSystem(sym_ham, t, q, p, sparams; simplify=false, cse=true))
    without = functions(HamiltonianSystem(sym_ham, t, q, p, sparams; simplify=false, cse=false))

    @test with.H(t₀, q₀, p₀, params) ≈ without.H(t₀, q₀, p₀, params) rtol = RTOL
    @test with.H(t₀, q₀, p₀, params) ≈ hamiltonian(t₀, q₀, p₀, params) rtol = RTOL

    for key in (:v, :f)
        compare(with, without, key, D, q₀, p₀)
    end
    compare(with, without, :ż, 2D, q₀, p₀)

    # `EHq = q̇ - ∂H/∂p` and `EHp = ṗ + ∂H/∂q` are residuals, so they take d/dt(q) and d/dt(p) in
    # trailing `V` and `F` slots.
    compare(with, without, :EHq, D, q₀, p₀, v₀, p₁)
    compare(with, without, :EHp, D, q₀, p₀, v₀, p₁)
    compare(with, without, :EH, 2D, q₀, p₀, v₀, p₁)
end
