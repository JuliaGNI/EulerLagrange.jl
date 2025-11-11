using EulerLagrange
using GeometricEquations
using LinearAlgebra
using Parameters
using Test


const tspan = (0.0, 3)
const tstep = 0.5

const G = 2.95912208286e-4

#Sun
const m₁ = 1.0
const q₁ = [0.0, 0.0, 0.0]
const q̇₁ = [0.0, 0.0, 0.0]
const p₁ = [0.0, 0.0, 0.0]

#Jupiter
const m₂ = 0.000954786104043
const q₂ = [-3.5023653, -3.8169847, -1.5507963]
const q̇₂ = [0.00565429, -0.0041249, -0.00190589]
const p₂ = m₂ * q̇₂

#Saturn
const m₃ = 0.000285583733151
const q₃ = [9.0755314, -3.0458353, -1.6483708]
const q̇₃ = [0.00168318, 0.00483525, 0.00192462]
const p₃ = m₃ * q̇₃

#Uranus
const m₄ = 0.0000437273164546
const q₄ = [8.3101420, -16.2901086, -7.2521278]
const q̇₄ = [0.00354178, 0.00137102, 0.00055029]
const p₄ = m₄ * q̇₄

#Neptune
const m₅ = 0.0000517759138449
const q₅ = [11.4707666, -25.7294829, -10.8169456]
const q̇₅ = [0.00288930, 0.00114527, 0.00039677]
const p₅ = m₅ * q̇₅

#Pluto
const m₆ = 1 / (1.3e8)
const q₆ = [-15.5387357, -25.2225594, -3.1902382]
const q̇₆ = [0.00276725, -0.00170702, -0.00136504]
const p₆ = m₆ * q̇₆

# 2d System
const ss2 = (
    d=3,
    n=2,
    m=[m₁, m₂],
    t₀=0.0,
    q₀=[q₁; q₂],
    v₀=[q̇₁; q̇₂],
    p₀=[p₁; p₂],
    default_parameters=(
        d=3,
        n=2,
        G=2.95912208286e-4,
        m=[m₁, m₂]
    )
)

# 6d System
const ss6 = (
    d=3,
    n=6,
    m=[m₁, m₂, m₃, m₄, m₅, m₆],
    t₀=0.0,
    q₀=[q₁; q₂; q₃; q₄; q₅; q₆],
    v₀=[q̇₁; q̇₂; q̇₃; q̇₄; q̇₅; q̇₆],
    p₀=[p₁; p₂; p₃; p₄; p₅; p₆],
    default_parameters=(
        d=3,
        n=6,
        G=2.95912208286e-4,
        m=[m₁, m₂, m₃, m₄, m₅, m₆]
    )
)


function v̄(v, t, q, p, params)
    p = transpose(reshape(p, params.d, params.n))
    v .= reshape(transpose(p ./ params.m), params.d * params.n)
    nothing
end

function p̃(p, t, q, v, params)
    v = transpose(reshape(v, params.d, params.n))
    p .= reshape(transpose(v .* params.m), params.d * params.n)
    nothing
end

function hamiltonian(t, q, p, params, d, n)
    @unpack G, m = params

    q = transpose(reshape(q, d, n))
    p = transpose(reshape(p, d, n))

    H = zero(eltype(q))
    for i in 1:n
        H += 1 / 2 / m[i] * p[i, :] ⋅ p[i, :]
        for j in 1:i-1
            H -= G * (m[i] * m[j]) / norm(q[i, :] - q[j, :])
        end
    end

    return H
end


function lagrangian(t, q, v, params, d, n)
    @unpack G, m = params

    q = transpose(reshape(q, d, n))
    v = transpose(reshape(v, d, n))

    L = zero(eltype(q))
    for i in 1:n
        L += 1 / 2 * m[i] * v[i, :] ⋅ v[i, :]
        for j in 1:i-1
            L += G * (m[i] * m[j]) / norm(q[i, :] - q[j, :])
        end
    end

    return L
end



function test_solar_system(ss)
    @unpack d, n, m, t₀, q₀, v₀, p₀ = ss
    params = ss.default_parameters

    # Symbolic variables and parameters
    t, x, v = lagrangian_variables(d * n)
    sparams = symbolize(params)

    # LagrangianSystem
    lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams, d, n), t, x, v, sparams; simplify=false)

    p₁, p₂ = zero(p₀), zero(p₀)
    ṗ₁, ṗ₂ = zero(p₀), zero(p₀)

    eqs = functions(lag_sys)

    eqs.ϑ(p₁, t₀, q₀, v₀, params)
    eqs.f(ṗ₁, t₀, q₀, v₀, params)

    p̃(p₂, t₀, q₀, v₀, params)
    # f̃(ṗ₂, t₀, q₀, v₀, params)

    @test eqs.L(t₀, q₀, v₀, params) ≈ lagrangian(t₀, q₀, v₀, params, d, n) atol = sqrt(2eps())
    @test p₁ ≈ p₂ atol = 2eps()
    # @test ṗ₁ ≈ ṗ₂  atol=2eps()

end


test_solar_system(ss2)
test_solar_system(ss6)
