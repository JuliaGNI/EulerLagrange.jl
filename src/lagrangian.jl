
function substitute_ẋ_with_v!(eqs, ẋ, v)
    for i in eachindex(eqs)
        eqs[i] = substitute(eqs[i], [ẋᵢ=>vᵢ for (ẋᵢ,vᵢ) in zip(ẋ,v)])
    end
end

function substitute_lagrangian_variables!(eqs, x, v, X, V)
    for i in eachindex(eqs)
        eqs[i] = substitute(eqs[i], [z=>Z for (z,Z) in zip([x..., v...], [X..., V...])])
    end
end


struct LagrangianSystem
    L
    t
    x
    v
    parameters
    equations

    function LagrangianSystem(L, t, x, v, params = NamedTuple())

        @assert eachindex(x) == eachindex(v)

        RuntimeGeneratedFunctions.init(@__MODULE__)

        @variables X[axes(x, 1)]
        @variables V[axes(v, 1)]
        @variables P[axes(v, 1)]
        @variables F[axes(x, 1)]

        Dt = Differential(t)
        Dx = collect(Differential.(x))
        Dv = collect(Differential.(v))
        Dz = vcat(Dx,Dv)
        ẋ  = collect(Dt.(x))

        EL = [expand_derivatives(Dx[i](L) - Dt(Dv[i](L))) for i in eachindex(Dx,Dv)]
        f  = [expand_derivatives(dx(L)) for dx in Dx]
        g  = [expand_derivatives(Dt(dv(L))) for dv in Dv]
        ϑ  = [expand_derivatives(dv(L)) for dv in Dv]
        θ  = vcat(ϑ, zero(ϑ))
        ω  = [expand_derivatives(simplify(Dz[i](θ[j]) - Dz[j](θ[i]))) for i in eachindex(Dz,θ), j in eachindex(Dz,θ)]
        Ω  = [expand_derivatives(simplify(Dx[i](ϑ[j]) - Dx[j](ϑ[i]))) for i in eachindex(Dx,ϑ), j in eachindex(Dx,ϑ)]
        M  = [expand_derivatives(simplify(Dv[i](ϑ[j]))) for i in eachindex(Dv), j in eachindex(ϑ)]
        N  = [expand_derivatives(simplify(Dx[i](ϑ[j]))) for i in eachindex(Dv), j in eachindex(ϑ)]

        σ  = simplify.(inv(ω))
        Σ  = simplify.(inv(Ω))
        
        for eq in (EL, f, g, ϑ, ω, Ω, M, N, σ, Σ)
            substitute_ẋ_with_v!(eq, ẋ, v)
            substitute_lagrangian_variables!(eq, x, v, X, V)
        end

        L = substitute(L, [z=>Z for (z,Z) in zip([x..., v...], [X..., V...])])
        a = inv(M) * (f - N * V)
        ϕ = P .- ϑ
        ψ = F .- g

        code_EL = build_function(EL, t, X, V)[2]
        code_a  = build_function(a,  t, X, V)[2]
        code_f  = build_function(f,  t, X, V)[2]
        code_g  = build_function(g,  t, X, V)[2]
        code_p  = build_function(ϑ,  t, X, V)[1]
        code_ϑ  = build_function(ϑ,  t, X, V)[2]
        code_ω  = build_function(ω,  t, X, V)[2]
        code_Ω  = build_function(Ω,  t, X, V)[2]
        code_ϕ  = build_function(ϕ,  t, X, V, P)[2]
        code_ψ  = build_function(ψ,  t, X, V, P, F)[2]
        code_L  = build_function(L,  t, X, V)
        code_M  = build_function(M,  t, X, V)[2]
        code_P  = build_function(Σ,  t, X, V)[2]

        eqs = (
            EL = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EL)),
            a  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_a)),
            f  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_f)),
            g  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_g)),
            p  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_p)),
            ϑ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ϑ)),
            ω  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ω)),
            Ω  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_Ω)),
            ϕ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ϕ)),
            ψ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ψ)),
            L  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_L)),
            M  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_M)),
            P  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_P)),
        )

        new(L, t, x, v, params, eqs)
    end
end

lagrangian(lsys::LagrangianSystem) = lsys.L
variables(lsys::LagrangianSystem) = (lsys.t, lsys.x, lsys.v)
parameters(lsys::LagrangianSystem) = lsys.parameters
equations(lsys::LagrangianSystem) = lsys.equations
equation(lsys::LagrangianSystem, name::Symbol) = equations(lsys)[name]


function lagrangian_variables(dimension::Int)
    t = parameter(:t)
    
    @variables (x(t))[1:dimension]
    @variables (v(t))[1:dimension]

    return (t,x,v)
end

function equation_wrappers(lsys::LagrangianSystem)
    (
        ϑ = (ϑ,t,x,v,params)   -> equation(lsys, :ϑ)(ϑ,t,x,v),
        f = (f,t,x,v,params)   -> equation(lsys, :f)(f,t,x,v),
        g = (f̄,t,x,v,λ,params) -> equation(lsys, :f̄)(f̄,t,x,λ),
        ω = (ω,t,x,v,params)   -> equation(lsys, :ω)(ω,t,x,v),
        l = (t,x,v,params)     -> equation(lsys, :L)(t,x,v),
        v̄ = (v,t,x,p,params)   -> equation(lsys, :ẋ)(v,t,x),
        f̄ = (f,t,x,v,params)   -> equation(lsys, :f)(f,t,x,v)
    )
end


function LODE(lsys::LagrangianSystem)
    eqs = equation_wrappers(lsys)
    LODE(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.l; v̄ = eqs.v̄, f̄ = eqs.f̄)
end

function LODEProblem(lsys::LagrangianSystem, tspan::Tuple, tstep::Real, ics::NamedTuple)
    eqs = equation_wrappers(lsys)
    LODEProblem(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.l, tspan, tstep, ics; v̄ = eqs.v̄, f̄ = eqs.f̄)
end

function LODEProblem(lsys::LagrangianSystem, tspan::Tuple, tstep::Real, q₀::State, p₀::State, λ₀::State = zero(q₀);)
    ics = (q = q₀, p = p₀, λ = λ₀)
    LODEProblem(lsys, tspan, tstep, ics)
end
