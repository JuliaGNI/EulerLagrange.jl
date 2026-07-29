
function substitute_hamiltonian_variables(equ, q, p)
    @variables Q[axes(q, 1)]
    @variables P[axes(p, 1)]
    substitute(equ, [zᵢ => Zᵢ for (zᵢ, Zᵢ) in zip([q..., p...], [Q..., P...])])
end

function substitute_hamiltonian_variables(equs::AbstractArray, q, p)
    [substitute_hamiltonian_variables(eq, q, p) for eq in equs]
end

function substitute_hamiltonian_variables(equs::NamedTuple, q, p)
    NamedTuple{keys(equs)}(Tuple(substitute_hamiltonian_variables(eq, q, p) for eq in equs))
end


"""
    substitute_time_derivatives(equ, dq, V, dp, F)

Replace the time derivatives `dq = d/dt(q)` and `dp = d/dt(p)` by the algebraic variables `V`, `F`.

The residual forms `EHq = q̇ - ∂H/∂p` and `EHp = ṗ + ∂H/∂q` contain those derivatives, and nothing
else in the pipeline removes them, so without this they reach `build_function` as derivatives of
variables the generated function never receives.
"""
function substitute_time_derivatives(equ, dq, V, dp, F)
    pairs = vcat([dqᵢ => Vᵢ for (dqᵢ, Vᵢ) in zip(dq, collect(V))],
        [dpᵢ => Fᵢ for (dpᵢ, Fᵢ) in zip(dp, collect(F))])
    substitute(equ, pairs)
end

function substitute_time_derivatives(equs::AbstractArray, dq, V, dp, F)
    [substitute_time_derivatives(eq, dq, V, dp, F) for eq in equs]
end

function substitute_time_derivatives(equs::NamedTuple, dq, V, dp, F)
    NamedTuple{keys(equs)}(Tuple(substitute_time_derivatives(eq, dq, V, dp, F) for eq in equs))
end

function hamiltonian_variables(dimension::Int)
    @variables t
    @variables q(t)[1:dimension]
    @variables p(t)[1:dimension]
    return (t, q, p)
end

function hamiltonian_derivatives(t, q, p)
    Dt = Differential(t)
    Dq = collect(Differential.(q))
    Dp = collect(Differential.(p))

    return (Dt, Dq, Dp)
end


"""
    HamiltonianSystem
"""
struct HamiltonianSystem
    H
    t
    q
    p
    parameters
    equations
    functions

    function HamiltonianSystem(H, t, q, p, params=NamedTuple(); simplify=true, scalarize=true, cse=true)

        @assert eachindex(q) == eachindex(p)

        @variables Q[axes(q, 1)]
        @variables P[axes(p, 1)]

        # Algebraic stand-ins for the time derivatives that appear in the residual form of the
        # equations of motion: `V` for d/dt(q) and `F` for d/dt(p).
        @variables V[axes(q, 1)]
        @variables F[axes(p, 1)]

        Dt, Dq, Dp = hamiltonian_derivatives(t, q, p)

        Hs = scalarize ? Symbolics.scalarize(H) : H
        Hs = simplify ? Symbolics.simplify(Hs) : Hs
        Hs = Num(Hs)

        EHq = [expand_derivatives(Dt(q[i]) - Dp[i](Hs)) for i in eachindex(Dp, q)]
        EHp = [expand_derivatives(Dt(p[i]) + Dq[i](Hs)) for i in eachindex(Dq, p)]
        EH = vcat(EHq, EHp)
        v = [expand_derivatives(dp(Hs)) for dp in Dp]
        f = [expand_derivatives(-dq(Hs)) for dq in Dq]
        ż = vcat(v, f)

        equs = (
            H=Hs,
            EH=EH,
            EHq=EHq,
            EHp=EHp,
            v=v,
            f=f,
            ż=ż,
        )

        # `EHq`, `EHp` and `EH` are residuals, so they carry the time derivatives d/dt(q) and d/dt(p)
        # explicitly. Replace those by the algebraic `V` and `F` before code generation: otherwise
        # they survive as derivatives of variables that do not exist in the generated function, which
        # makes those three throw as soon as they are called.
        equs = substitute_time_derivatives(equs, collect(Dt.(q)), V, collect(Dt.(p)), F)

        equs_subs = substitute_hamiltonian_variables(equs, q, p)

        _build(expr, extra...) = build_function(expr, t, Q, P, extra..., params...; nanmath=false, cse=cse)

        code = (
            H=substitute_parameters(_build(equs_subs.H), params),
            EH=substitute_parameters(_build(equs_subs.EH, V, F)[2], params),
            EHq=substitute_parameters(_build(equs_subs.EHq, V, F)[2], params),
            EHp=substitute_parameters(_build(equs_subs.EHp, V, F)[2], params),
            v=substitute_parameters(_build(equs_subs.v)[2], params),
            f=substitute_parameters(_build(equs_subs.f)[2], params),
            ż=substitute_parameters(_build(equs_subs.ż)[2], params),
        )

        funcs = generate_code(code)

        new(Hs, t, q, p, params, equs, funcs)
    end
end

hamiltonian(hsys::HamiltonianSystem) = hsys.H
parameters(hsys::HamiltonianSystem) = hsys.parameters
variables(hsys::HamiltonianSystem) = (hsys.t, hsys.q, hsys.p)
equations(lsys::HamiltonianSystem) = lsys.equations
functions(lsys::HamiltonianSystem) = lsys.functions

function Base.show(io::IO, hsys::HamiltonianSystem)
    print(io, "\nHamiltonian system with\n")
    print(io, "\nH = ")
    print(io, hamiltonian(hsys))
    # print(io, "\n\nand equations of motion\n\n")
    # for eq in equations(hsys).EH
    #     print(io, eq)
    #     print(io, "\n")
    # end
end


function HODE(lsys::HamiltonianSystem; kwargs...)
    eqs = functions(lsys)
    HODE(eqs.v, eqs.f, eqs.H; kwargs...)
end

function HODEProblem(lsys::HamiltonianSystem, tspan::Tuple, tstep::Real, ics...; kwargs...)
    eqs = functions(lsys)
    HODEProblem(eqs.v, eqs.f, eqs.H, tspan, tstep, ics...; kwargs...)
end
