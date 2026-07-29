using EulerLagrange
using LinearAlgebra
using Test

# The generated functions are `RuntimeGeneratedFunction`s, which silently ignore trailing arguments.
# A signature that disagrees with its caller in GeometricEquations therefore does not throw: it reads
# the wrong values out of the wrong slots and drops the rest. `check_methods` does not catch it
# either, since `applicable` is satisfied at any arity. So pin the emitted signatures directly.
#
# Each entry below records the argument names `build_function` emitted, `ˍ₋out` being the in-place
# output buffer and `params` the parameter tuple appended by `substitute_parameters`. The callers
# these must satisfy live in GeometricEquations: `_get_*(::LODE, params)` in `src/odes/lode.jl` and
# `_get_*(::LDAE, params)` in `src/daes/ldae.jl`.

argnames(f) = typeof(f).parameters[1]

function test_signatures(funcs, expected)
    @test Set(keys(funcs)) == Set(keys(expected))
    for (key, names) in pairs(expected)
        @test argnames(getproperty(funcs, key)) == names
    end
end

t, x, v = lagrangian_variables(2)


@testset "$(rpad("LagrangianSystem", 80))" begin
    L = (v ⋅ v) / 2 - (x ⋅ x) / 2

    test_signatures(functions(LagrangianSystem(L, t, x, v)), (
        L=(:t, :X, :V, :params),
        # `Λ` is the acceleration d/dt(v). `g` is called by GeometricEquations as
        # `g(out, t, q, v, λ, params)`, so `V` must precede `Λ`.
        EL=(:ˍ₋out, :t, :X, :V, :Λ, :params),
        f=(:ˍ₋out, :t, :X, :V, :params),
        g=(:ˍ₋out, :t, :X, :V, :Λ, :params),
        p=(:t, :X, :V, :params),
        ϑ=(:ˍ₋out, :t, :X, :V, :params),
        # 2n×2n for a regular Lagrangian; called as `ω(out, t, q, v, params)`.
        ω=(:ˍ₋out, :t, :X, :V, :params),
        ϕ=(:ˍ₋out, :t, :X, :V, :P, :params),
        ψ=(:ˍ₋out, :t, :X, :V, :P, :F, :Λ, :params),
        M=(:ˍ₋out, :t, :X, :V, :params),
    ))
end


@testset "$(rpad("DegenerateLagrangianSystem", 80))" begin
    K = (log(x[2]) / x[1] / 2) * v[1] - (log(x[1]) / x[2] / 2) * v[2]
    H = x[1] + x[2] + log(x[1]) + log(x[2])

    test_signatures(functions(DegenerateLagrangianSystem(K, H, t, x, v)), (
        L=(:t, :X, :V, :params),
        H=(:t, :X, :params),
        EL=(:ˍ₋out, :t, :X, :V, :params),
        ∇H=(:ˍ₋out, :t, :X, :V, :params),
        ẋ=(:ˍ₋out, :t, :X, :params),
        v=(:ˍ₋out, :t, :X, :P, :params),
        f=(:ˍ₋out, :t, :X, :V, :params),
        # `u`, `g`, `ū` and `ḡ` are reached through the wrappers in `lagrangian_degenerate.jl`,
        # which map LDAE's `(t, q, v, p, λ)` onto these slots: the multiplier λ arrives in `V`.
        u=(:ˍ₋out, :t, :X, :Λ, :V, :params),
        g=(:ˍ₋out, :t, :X, :Λ, :V, :params),
        ū=(:ˍ₋out, :t, :X, :Λ, :P, :V, :params),
        ḡ=(:ˍ₋out, :t, :X, :Λ, :P, :V, :params),
        p=(:t, :X, :V, :params),
        ϑ=(:ˍ₋out, :t, :X, :V, :params),
        # n×n for a degenerate Lagrangian, which is already a first-order system.
        ω=(:ˍ₋out, :t, :X, :V, :params),
        ϕ=(:ˍ₋out, :t, :X, :V, :P, :params),
        # `_get_ψ(::LDAE, params)` calls `ψ(out, t, q, v, p, q̇, ṗ, params)`: `Ẋ` is q̇ and `F` is ṗ.
        ψ=(:ˍ₋out, :t, :X, :V, :P, :Ẋ, :F, :params),
        P=(:ˍ₋out, :t, :X, :V, :params),
    ))
end


@testset "$(rpad("HamiltonianSystem", 80))" begin
    tₕ, q, p = hamiltonian_variables(2)
    H = (p ⋅ p) / 2 + (q ⋅ q) / 2

    # `V` and `F` stand in for d/dt(q) and d/dt(p) in the residual forms.
    test_signatures(functions(HamiltonianSystem(H, tₕ, q, p)), (
        H=(:t, :Q, :P, :params),
        EH=(:ˍ₋out, :t, :Q, :P, :V, :F, :params),
        EHq=(:ˍ₋out, :t, :Q, :P, :V, :F, :params),
        EHp=(:ˍ₋out, :t, :Q, :P, :V, :F, :params),
        v=(:ˍ₋out, :t, :Q, :P, :params),
        f=(:ˍ₋out, :t, :Q, :P, :params),
        ż=(:ˍ₋out, :t, :Q, :P, :params),
    ))
end
