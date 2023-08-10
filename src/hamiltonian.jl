function substitute_hamiltonian_variables!(eqs, q, p, Q, P)
    for i in eachindex(eqs)
        eqs[i] = substitute(eqs[i], [z=>Z for (z,Z) in zip([q..., p...], [Q..., P...])])
    end
end

struct HamiltonianSystem
    t
    q
    p
    H
    # parameters
    # functions
    eqs

    function HamiltonianSystem(t, q, p, H)

        @assert eachindex(q) == eachindex(p)

        RuntimeGeneratedFunctions.init(@__MODULE__)

        @variables Q[axes(q,1)]
        @variables P[axes(p,1)]

        Dt = Differential(t)
        Dq = collect(Differential.(q))
        Dp = collect(Differential.(p))

        EHq = [expand_derivatives(Dq[i](H) - Dt(p[i])) for i in eachindex(Dq,p)]
        EHp = [expand_derivatives(Dp[i](H) + Dt(q[i])) for i in eachindex(Dp,q)]
        EH  = vcat(EHq, EHp)
        f   = [expand_derivatives(-dq(H)) for dq in Dq]
        v   = [expand_derivatives( dp(H)) for dp in Dp]
        ż   = vcat(-v, -f)

        for eq in (EH, EHq, EHp, f, v, ż)
            substitute_hamiltonian_variables!(eq, q, p, Q, P)
        end

        H = substitute(H, [z=>Z for (z,Z) in zip([q..., p...], [Q..., P...])])

        code_EH  = build_function(EH, t, Q, P)[1]
        code_EHq = build_function(EHq, t, Q, P)[1]
        code_EHp = build_function(EHp, t, Q, P)[2]
        code_H   = build_function(H, t, Q, P)
        code_f   = build_function(f, t, Q, P)[1]
        code_v   = build_function(v, t, Q, P)[1]
        code_ż   = build_function(ż, t, Q, P)[1]

        eqs = (
            EH  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EH)),
            EHq = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EHq)),
            EHp = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EHp)),
            H   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_H)),
            f   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_f)),
            v   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_v)),
            ż   = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ż))
        )

        new(t, q, p, H, eqs)

    end

end


function hamiltonian_variables(dimension::Int)
    t = parameter(:t)
    @variables (q(t))[1:dimension]
    @variables (p(t))[1:dimension]
    return (t,q,p)
end
