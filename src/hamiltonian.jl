
function substitute_hamiltonian_variables(equ, q, p, Q, P)
    substitute(equ, [z=>Z for (z,Z) in zip([q..., p...], [Q..., P...])])
end

function substitute_hamiltonian_variables!(eqs, q, p, Q, P)
    for i in eachindex(eqs)
        eqs[i] = substitute_hamiltonian_variables(eqs[i], q, p, Q, P)
    end
end

struct HamiltonianSystem
    H
    t
    q
    p
    parameters
    equations

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
            substitute_hamiltonian_variables!(eq, q, p, Q, P)
        end

        H  = substitute_hamiltonian_variables(H, q, p, Q, P)
        EH = vcat(EHq, EHp)
        ż  = vcat(v, f)

        code_EH  = substitute_parameters(build_function(EH, t, Q, P, params...)[2],  params)
        code_EHq = substitute_parameters(build_function(EHq, t, Q, P, params...)[2], params)
        code_EHp = substitute_parameters(build_function(EHp, t, Q, P, params...)[2], params)
        code_H   = substitute_parameters(build_function(H, t, Q, P, params...),      params)
        code_v   = substitute_parameters(build_function(v, t, Q, P, params...)[2],   params)
        code_f   = substitute_parameters(build_function(f, t, Q, P, params...)[2],   params)
        code_ż   = substitute_parameters(build_function(ż, t, Q, P, params...)[2],   params)

        eqs = (
            EH  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EH)),
            EHq = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EHq)),
            EHp = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EHp)),
            H   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_H)),
            v   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_v)),
            f   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_f)),
            ż   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ż))
        )

        param_eqs = length(params) > 0 ? eqs : (
            EH  = (EH, t,q,p,params) -> eqs.EH(EH,t,q,p),
            EHq = (EHq,t,q,p,params) -> eqs.EHq(EHq,t,q,p),
            EHp = (EHp,t,q,p,params) -> eqs.EHp(EHp,t,q,p),
            H = (t,q,p,params)       -> eqs.H(t,q,p),
            v = (v,t,q,p,params)     -> eqs.v(v,t,q,p),
            f = (f,t,q,p,params)     -> eqs.f(f,t,q,p),
            ż = (ż,t,q,p,params)     -> eqs.ż(ż,t,q,p),
        )

        new(H, t, q, p, params, param_eqs)
    end
end

hamiltonian(hsys::HamiltonianSystem) = hsys.H
variables(hsys::HamiltonianSystem) = (hsys.t, hsys.q, hsys.p)
parameters(hsys::HamiltonianSystem) = hsys.parameters
equations(hsys::HamiltonianSystem) = hsys.equations
equation(hsys::HamiltonianSystem, name::Symbol) = equations(hsys)[name]


function hamiltonian_variables(dimension::Int)
    t = parameter(:t)
    @variables (q(t))[1:dimension]
    @variables (p(t))[1:dimension]
    return (t,q,p)
end
