
function substitute_hamiltonian_variables(equ, q, p)
    @variables Q[axes(q,1)]
    @variables P[axes(p,1)]
    substitute(equ, [zᵢ=>Zᵢ for (zᵢ,Zᵢ) in zip([q..., p...], [Q..., P...])])
end

function substitute_hamiltonian_variables(equs::Union{AbstractArray, ArrayLike}, q, p)
    [substitute_hamiltonian_variables(eq, q, p) for eq in equs]
end

function substitute_hamiltonian_variables(equs::NamedTuple, q, p)
    NamedTuple{keys(equs)}(Tuple(substitute_hamiltonian_variables(eq, q, p) for eq in equs))
end

struct HamiltonianSystem
    H
    t
    q
    p
    parameters
    equations
    functions

    function HamiltonianSystem(H, t, q, p, params = NamedTuple())

        @assert eachindex(q) == eachindex(p)

        @variables Q[axes(q,1)]
        @variables P[axes(p,1)]

        Dt = Differential(t)
        Dq = collect(Differential.(q))
        Dp = collect(Differential.(p))

        EHq = [expand_derivatives(Dt(q[i]) - Dp[i](H)) for i in eachindex(Dp,q)]
        EHp = [expand_derivatives(Dt(p[i]) + Dq[i](H)) for i in eachindex(Dq,p)]
        EH  = vcat(EHq, EHp)
        v   = [expand_derivatives( dp(H)) for dp in Dp]
        f   = [expand_derivatives(-dq(H)) for dq in Dq]
        ż   = vcat(v, f)

        equs = (
            H = H,
            EH = EH,
            EHq = EHq,
            EHp = EHp,
            v = v,
            f = f,
            ż = ż,
        )

        equs_subs = substitute_hamiltonian_variables(equs, q, p)

        code = (
            H   = substitute_parameters(build_function(equs_subs.H, t, Q, P, params...),      params),
            EH  = substitute_parameters(build_function(equs_subs.EH, t, Q, P, params...)[2],  params),
            EHq = substitute_parameters(build_function(equs_subs.EHq, t, Q, P, params...)[2], params),
            EHp = substitute_parameters(build_function(equs_subs.EHp, t, Q, P, params...)[2], params),
            v   = substitute_parameters(build_function(equs_subs.v, t, Q, P, params...)[2],   params),
            f   = substitute_parameters(build_function(equs_subs.f, t, Q, P, params...)[2],   params),
            ż   = substitute_parameters(build_function(equs_subs.ż, t, Q, P, params...)[2],   params),
        )

        funcs = generate_code(code)

        new(H, t, q, p, params, equs, funcs)
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
    print(io, "\n\nand equations of motion\n\n")
    for eq in equations(hsys).EH
        print(io, eq)
        print(io, "\n")
    end
end

function hamiltonian_variables(dimension::Int)
    t = parameter(:t)
    @variables (q(t))[1:dimension]
    @variables (p(t))[1:dimension]
    return (t,q,p)
end


function HODE(lsys::HamiltonianSystem; kwargs...)
    eqs = functions(lsys)
    HODE(eqs.v, eqs.f, eqs.H; kwargs...)
end

function HODEProblem(lsys::HamiltonianSystem, tspan::Tuple, tstep::Real, ics::NamedTuple; kwargs...)
    eqs = functions(lsys)
    HODEProblem(eqs.v, eqs.f, eqs.H, tspan, tstep, ics; kwargs...)
end

function HODEProblem(lsys::HamiltonianSystem, tspan::Tuple, tstep::Real, q₀::StateVariable, p₀::StateVariable; kwargs...)
    ics = (q = q₀, p = p₀)
    HODEProblem(lsys, tspan, tstep, ics; kwargs...)
end
