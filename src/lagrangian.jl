
function substitute_ẋ_with_v(equ, ẋ, v)
    substitute(equ, [ẋᵢ=>vᵢ for (ẋᵢ,vᵢ) in zip(ẋ,v)])
end

function substitute_ẋ_with_v!(eqs, ẋ, v)
    for i in eachindex(eqs)
        eqs[i] = substitute_ẋ_with_v(eqs[i], ẋ, v)
    end
end

function substitute_lagrangian_variables(equ, x, v, X, V)
    substitute(equ, [z=>Z for (z,Z) in zip([x..., v...], [X..., V...])])
end

function substitute_lagrangian_variables!(eqs, x, v, X, V)
    for i in eachindex(eqs)
        eqs[i] = substitute_lagrangian_variables(eqs[i], x, v, X, V)
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

        L = substitute_lagrangian_variables(L, x, v, X, V)
        a = inv(M) * (f - N * V)
        ϕ = P .- ϑ
        ψ = F .- g

        code_EL = substitute_parameters(build_function(EL, t, X, V, params...)[2], params)
        code_a  = substitute_parameters(build_function(a,  t, X, V, params...)[2], params)
        code_f  = substitute_parameters(build_function(f,  t, X, V, params...)[2], params)
        code_g  = substitute_parameters(build_function(g,  t, X, V, params...)[2], params)
        code_p  = substitute_parameters(build_function(ϑ,  t, X, V, params...)[1], params)
        code_ϑ  = substitute_parameters(build_function(ϑ,  t, X, V, params...)[2], params)
        code_ω  = substitute_parameters(build_function(ω,  t, X, V, params...)[2], params)
        code_Ω  = substitute_parameters(build_function(Ω,  t, X, V, params...)[2], params)
        code_ϕ  = substitute_parameters(build_function(ϕ,  t, X, V, P, params...)[2], params)
        code_ψ  = substitute_parameters(build_function(ψ,  t, X, V, P, F, params...)[2], params)
        code_L  = substitute_parameters(build_function(L,  t, X, V, params...), params)
        code_M  = substitute_parameters(build_function(M,  t, X, V, params...)[2], params)
        code_P  = substitute_parameters(build_function(Σ,  t, X, V, params...)[2], params)

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

        param_eqs = length(params) > 0 ? eqs : (
            L = (t,x,v,params)     -> eqs.L(t,x,v),
            EL= (e,t,x,v,params)   -> eqs.EL(e,t,x,v),
            a = (a,t,x,v,params)   -> eqs.a(a,t,x,v),
            ϑ = (ϑ,t,x,v,params)   -> eqs.ϑ(ϑ,t,x,v),
            f = (f,t,x,v,params)   -> eqs.f(f,t,x,v),
            g = (f̄,t,x,v,λ,params) -> eqs.f̄(f̄,t,x,λ),
            ω = (ω,t,x,v,params)   -> eqs.ω(ω,t,x,v),
            Ω = (Ω,t,x,v,params)   -> eqs.Ω(Ω,t,x,v),
            ϕ = (ϕ,t,x,v,params)   -> eqs.ϕ(ϕ,t,x,v),
            ψ = (ψ,t,x,v,params)   -> eqs.ψ(ψ,t,x,v),
            v̄ = (v,t,x,p,params)   -> eqs.ẋ(v,t,x,p),
            f̄ = (f,t,x,v,params)   -> eqs.f(f,t,x,v),
            M = (M,t,x,v,params)   -> eqs.M(M,t,x,v),
            P = (P,t,x,v,params)   -> eqs.P(P,t,x,v),
        )

        new(L, t, x, v, params, param_eqs)
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


function LODE(lsys::LagrangianSystem)
    eqs = equations(lsys)
    LODE(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L; v̄ = eqs.v̄, f̄ = eqs.f̄)
end

function LODEProblem(lsys::LagrangianSystem, tspan::Tuple, tstep::Real, ics::NamedTuple)
    eqs = equations(lsys)
    LODEProblem(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, tspan, tstep, ics; v̄ = eqs.v̄, f̄ = eqs.f̄)
end

function LODEProblem(lsys::LagrangianSystem, tspan::Tuple, tstep::Real, q₀::State, p₀::State, λ₀::State = zero(q₀);)
    ics = (q = q₀, p = p₀, λ = λ₀)
    LODEProblem(lsys, tspan, tstep, ics)
end
