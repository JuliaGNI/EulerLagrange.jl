@doc raw"""
    LagrangianSystem(L, t, x, v, params = NamedTuple(); simplify = false, scalarize = true,
                     cse = true, nanmath = false)

The equations of motion of a regular Lagrangian `L`, generated symbolically.

Since `L` is regular, the system is second order in ``n`` equations, equivalently first order in
``2n``; accordingly `ω` is the ``2n × 2n`` two-form on ``(q, \dot q)``. See
[`DegenerateLagrangianSystem`](@ref) for the first-order, ``n × n`` case.

# Keyword arguments

  - `simplify`: apply `Symbolics.simplify` to `L` before differentiating. **Off by default**, because
    it is at best neutral and often much worse: it rewrites a sum of fractions into a
    common-denominator form whose derivatives are far more expensive both to build and to evaluate.
    Measured across every EulerLagrange-based problem in GeometricProblems (88 generated functions),
    `simplify = true` was never faster to evaluate and up to 15× slower, at 23× the total
    construction cost.
  - `scalarize`: apply `Symbolics.scalarize` to `L`, expanding array expressions into components.
  - `cse`: eliminate common subexpressions in the generated code. On by default: measured over 45
    generated functions in GeometricProblems it was faster to evaluate 18 times, slower once, and
    equal otherwise, with the wins concentrated in the forces (up to 8.4× on the 18-degree-of-freedom
    N-body force). Note that it emits *larger* code in most cases — it binds constant nodes to
    temporaries as well as shared subexpressions — so size is not a proxy for speed here. Binding
    constants also means results may differ from `cse = false` in the last bit.
  - `nanmath`: emit `NaNMath` variants of functions like `log`, `sqrt` and `^`, which return `NaN`
    outside their real domain instead of throwing a `DomainError`. Off by default, so the generated
    code uses the ordinary `Base` functions and an out-of-domain state is reported as an error rather
    than propagating silently.
"""
struct LagrangianSystem
    L
    t
    x
    v
    parameters
    equations
    functions

    function LagrangianSystem(L, t, x, v, params=NamedTuple(); simplify=false, scalarize=true, cse=true, nanmath=false)

        @assert eachindex(x) == eachindex(v)

        @variables p(t)[axes(x, 1)]
        @variables f(t)[axes(x, 1)]

        @variables X[axes(x, 1)]
        @variables V[axes(v, 1)]
        @variables P[axes(v, 1)]
        @variables F[axes(x, 1)]
        @variables Λ[axes(v, 1)]

        Dt, Dx, Dv = lagrangian_derivatives(t, x, v)

        Dz = vcat(Dx, Dv)
        ẋ = collect(Dt.(x))
        ṗ = collect(Dt.(p))

        Ls = scalarize ? Symbolics.scalarize(L) : L
        Ls = simplify ? Symbolics.simplify(Ls) : Ls
        Ls = Num(Ls)

        f = [expand_derivatives(dx(Ls)) for dx in Dx]
        g = [expand_derivatives(Dt(dv(Ls))) for dv in Dv]
        ϑ = [expand_derivatives(dv(Ls)) for dv in Dv]
        EL = [f[i] - g[i] for i in eachindex(f, g)]

        # The Lagrange one-form θ_L = ϑᵢ dqⁱ expressed in the (x, v) coordinates that `Dz` spans: it
        # has the components of ϑ along dx and none along dv.
        #
        # `ω = dθ_L` is built from `ϑ = ∂L/∂v`, not from the full gradient `∂L/∂z`. Taking the latter
        # gives `ω = d(dL) ≡ 0` — an identically vanishing matrix, and formerly by far the most
        # expensive quantity in this constructor.
        #
        # A `LagrangianSystem` describes a *regular* Lagrangian: a second-order system of n equations,
        # equivalent to a first-order system of 2n. Its `ω` is therefore the full 2n×2n two-form on
        # (q, q̇), in block form `[N - Nᵀ  -Mᵀ; M  0]`, which reduces to the canonical `[0 -M; M 0]`
        # when ϑ depends on the velocities alone. A `DegenerateLagrangianSystem` is a first-order
        # system of n equations and its `ω` is correspondingly n×n; see `lagrangian_degenerate.jl`.
        ϑz = vcat(ϑ, zero(ϑ))

        ω = [Dz[i](ϑz[j]) - Dz[j](ϑz[i]) for i in eachindex(Dz, ϑz), j in eachindex(Dz, ϑz)]
        M = [Dv[i](ϑ[j]) for i in eachindex(Dv), j in eachindex(ϑ)]
        N = [Dx[i](ϑ[j]) for i in eachindex(Dv), j in eachindex(ϑ)]

        # `simplify` deliberately does not run here. At this point the entries are still unevaluated
        # `Differential` applications, so no rewrite rule can cancel anything — simplifying first
        # costs a great deal and leaves the expanded result unsimplified anyway. Simplifying *after*
        # the expansion is intractable at the sizes that matter, so `cse=true` in `build_function`
        # below is what keeps the emitted code small instead.
        ω = expand_derivatives.(ω)
        M = expand_derivatives.(M)
        N = expand_derivatives.(N)

        equs = (
            L=Ls,
            EL=EL,
            f=f,
            g=g,
            ϑ=ϑ,
            ω=ω,
            M=M,
            N=N,
        )

        # _simplify(expr, dosimplify) = dosimplify ? simplify.(expr) : expr

        # `g` — and so `EL = f - g` and `ψ = ṗ - g` — contains the acceleration d/dt(v). Replace it
        # by `Λ` before the remaining substitutions, otherwise it survives into the generated code as
        # a derivative of a variable that does not exist there. See `substitute_acceleration`.
        equs = substitute_acceleration(equs, collect(Dt.(v)), Λ)

        equs_subs = substitute_lagrangian_variables(equs, x, ẋ, v)
        equs_subs = merge(equs_subs, (
            # a = inv(equs_subs.M) * (equs_subs.f - equs_subs.N * V),
            ϕ=P .- equs_subs.ϑ,
            ψ=F .- equs_subs.g,
            # σ = _simplify(inv(equs_subs.ω), dosimplify),
            # Σ = _simplify(inv(equs_subs.Ω), dosimplify),
        ))

        equs = substitute_v_with_ẋ(equs, v, ẋ)
        equs = merge(equs, (
            ϕ=p .- equs.ϑ,
            ψ=ṗ .- equs.g,
        ))

        _build(expr, args...) = build_function(expr, args...; nanmath=nanmath, cse=cse)

        # `p` and `ϑ` are the out-of-place and in-place halves of the same pair, so build once.
        ϑcode = _build(equs_subs.ϑ, t, X, V, params...)

        code = (
            L=substitute_parameters(_build(equs_subs.L, t, X, V, params...), params),
            # `EL`, `g` and `ψ` take the acceleration `Λ`. `g`'s argument order is fixed by
            # GeometricEquations, whose `_get_g(::LODE, params)` calls `g(out, t, q, v, λ, params)`;
            # it used to be built as `(t, X, Λ, V, …)`, with `Λ` and `V` the wrong way round.
            EL=substitute_parameters(_build(equs_subs.EL, t, X, V, Λ, params...)[2], params),
            # a  = substitute_parameters(_build(equs_subs.a,  t, X, V, params...)[2], params),
            f=substitute_parameters(_build(equs_subs.f, t, X, V, params...)[2], params),
            g=substitute_parameters(_build(equs_subs.g, t, X, V, Λ, params...)[2], params),
            p=substitute_parameters(ϑcode[1], params),
            ϑ=substitute_parameters(ϑcode[2], params),
            ω=substitute_parameters(_build(equs_subs.ω, t, X, V, params...)[2], params),
            ϕ=substitute_parameters(_build(equs_subs.ϕ, t, X, V, P, params...)[2], params),
            ψ=substitute_parameters(_build(equs_subs.ψ, t, X, V, P, F, Λ, params...)[2], params),
            M=substitute_parameters(_build(equs_subs.M, t, X, V, params...)[2], params),
            # P  = substitute_parameters(_build(equs_subs.Σ,  t, X, V, params...)[2], params),
        )

        new(Ls, t, x, v, params, equs, generate_code(code))
    end
end

lagrangian(lsys::LagrangianSystem) = lsys.L
parameters(lsys::LagrangianSystem) = lsys.parameters
variables(lsys::LagrangianSystem) = (lsys.t, lsys.x, lsys.v)
equations(lsys::LagrangianSystem) = lsys.equations
functions(lsys::LagrangianSystem) = lsys.functions

function Base.show(io::IO, lsys::LagrangianSystem)
    print(io, "\nLagrangian system with\n")
    print(io, "\nL = ")
    print(io, lagrangian(lsys))
    # print(io, "\n\nand equations of motion\n\n")
    # for eq in equations(lsys).EL
    #     print(io, eq)
    #     print(io, "\n")
    # end
end


function LODE(lsys::LagrangianSystem; kwargs...)
    eqs = functions(lsys)
    LODE(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L; kwargs...)
end

function LODEProblem(lsys::LagrangianSystem, tspan::Tuple, tstep::Real, ics...; kwargs...)
    eqs = functions(lsys)
    LODEProblem(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, tspan, tstep, ics...; kwargs...)
end
