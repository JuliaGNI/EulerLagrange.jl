
function substitute_hamiltonian_variables(equ, q, p)
    @variables Q[axes(q,1)]
    @variables P[axes(p,1)]
    substitute(equ, [z=>Z for (z,Z) in zip([q..., p...], [Q..., P...])])
end

function substitute_hamiltonian_variables!(eqs, q, p)
    for i in eachindex(eqs)
        eqs[i] = substitute_hamiltonian_variables(eqs[i], q, p)
    end
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

        RuntimeGeneratedFunctions.init(@__MODULE__)

        @variables Q[axes(q,1)]
        @variables P[axes(p,1)]

        Dt = Differential(t)
        Dq = collect(Differential.(q))
        Dp = collect(Differential.(p))

        EHq = [expand_derivatives(Dt(q[i]) - Dp[i](H)) for i in eachindex(Dp,q)]
        EHp = [expand_derivatives(Dt(p[i]) + Dq[i](H)) for i in eachindex(Dq,p)]
        v   = [expand_derivatives( dp(H)) for dp in Dp]
        f   = [expand_derivatives(-dq(H)) for dq in Dq]

        for eq in (EHq, EHp, v, f)
            substitute_hamiltonian_variables!(eq, q, p)
        end

        H  = substitute_hamiltonian_variables(H, q, p)
        EH = vcat(EHq, EHp)
        ż  = vcat(v, f)

        equs = (
            H = H,
            EH = EH,
            EHq = EHq,
            EHp = EHp,
            v = v,
            f = f,
            ż = ż,
        )

        code = (
            H   = substitute_parameters(build_function(equs.H, t, Q, P, params...),      params),
            EH  = substitute_parameters(build_function(equs.EH, t, Q, P, params...)[2],  params),
            EHq = substitute_parameters(build_function(equs.EHq, t, Q, P, params...)[2], params),
            EHp = substitute_parameters(build_function(equs.EHp, t, Q, P, params...)[2], params),
            v   = substitute_parameters(build_function(equs.v, t, Q, P, params...)[2],   params),
            f   = substitute_parameters(build_function(equs.f, t, Q, P, params...)[2],   params),
            ż   = substitute_parameters(build_function(equs.ż, t, Q, P, params...)[2],   params),
        )

        funcs = NamedTuple{keys(code)}(Tuple(@RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(c)) for c in code))

        # funcs = (
        #     H   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_H)),
        #     EH  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EH)),
        #     EHq = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EHq)),
        #     EHp = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EHp)),
        #     v   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_v)),
        #     f   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_f)),
        #     ż   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ż))
        # )

        funcs_param = length(params) > 0 ? funcs : (
            H = (t,q,p,params)       -> funcs.H(t,q,p),
            EH  = (EH, t,q,p,params) -> funcs.EH(EH,t,q,p),
            EHq = (EHq,t,q,p,params) -> funcs.EHq(EHq,t,q,p),
            EHp = (EHp,t,q,p,params) -> funcs.EHp(EHp,t,q,p),
            v = (v,t,q,p,params)     -> funcs.v(v,t,q,p),
            f = (f,t,q,p,params)     -> funcs.f(f,t,q,p),
            ż = (ż,t,q,p,params)     -> funcs.ż(ż,t,q,p),
        )

        new(H, t, q, p, params, equs, funcs_param)
    end
end

hamiltonian(hsys::HamiltonianSystem) = hsys.H
parameters(hsys::HamiltonianSystem) = hsys.parameters
variables(hsys::HamiltonianSystem) = (hsys.t, hsys.q, hsys.p)
equations(lsys::HamiltonianSystem) = lsys.equations
functions(lsys::HamiltonianSystem) = lsys.functions


function hamiltonian_variables(dimension::Int)
    t = parameter(:t)
    @variables (q(t))[1:dimension]
    @variables (p(t))[1:dimension]
    return (t,q,p)
end


function HODE(lsys::HamiltonianSystem)
    eqs = functions(lsys)
    HODE(eqs.v, eqs.f, eqs.H)
end

function HODEProblem(lsys::HamiltonianSystem, tspan::Tuple, tstep::Real, ics::NamedTuple)
    eqs = functions(lsys)
    HODEProblem(eqs.v, eqs.f, eqs.H, tspan, tstep, ics)
end

function HODEProblem(lsys::HamiltonianSystem, tspan::Tuple, tstep::Real, q₀::State, p₀::State)
    ics = (q = q₀, p = p₀)
    HODEProblem(lsys, tspan, tstep, ics)
end
