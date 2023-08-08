function substitute_hamiltonian_variables!(eqs, q, p, Q, P)
    for i in eachindex(eqs)
        eqs[i] = substitute(eqs[i], [z=>Z for (z,Z) in zip([q..., p...], [Q..., P...])])
    end
end

struct HamiltonianSystem
    t
    q
    p
    params
    H
    # parameters
    # functions
    eqs

    function HamiltonianSystem(t, q, p, params, H)

        @assert eachindex(q) == eachindex(p)

        RuntimeGeneratedFunctions.init(@__MODULE__)

        @variables Q[eachindex(q)]
        @variables P[eachindex(p)]

        Dt = Differential(t)
        Dq = Differential.(q)
        Dp = Differential.(p)

        EHq = [expand_derivatives(Dq[i](H) - Dt(p[i])) for i in eachindex(Dq,p)]
        EHp = [expand_derivatives(Dp[i](H) + Dt(q[i])) for i in eachindex(Dp,q)]
        EH  = vcat(EHq, EHp)

        for eq in (EH, EHq, EHp)
            substitute_hamiltonian_variables!(eq, q, p, Q, P)
        end

        H = substitute(L, [z=>Z for (z,Z) in zip([q..., p...], [Q..., P...])])

        code_EH  = build_function(EH, t, Q, P)
        code_EHq = build_function(EHq, t, Q, P)[2]
        code_EHp = build_function(EHp, t, Q, P)[2]
        code_H   = build_function(H, t, Q, P)

        eqs = (
            EH  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EH)),
            EHq = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EHq)),
            EHp = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EHp)),
            H   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_H)),
        )

        new(t, q, p, params, H, eqs)

    end

end




function hamiltonian_variables(dimension::Int)
    t = parameter(:t)
    q = @variables (q(t))[1:dimension]
    p = @variables (p(t))[1:dimension]
    return (t,q,p)
end
