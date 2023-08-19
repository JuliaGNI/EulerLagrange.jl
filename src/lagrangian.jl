
function substitute_ẋ_with_v(equ, ẋ, v)
    substitute(equ, [ẋᵢ=>vᵢ for (ẋᵢ,vᵢ) in zip(ẋ,v)])
end

function substitute_ẋ_with_v!(eqs, ẋ, v)
    for i in eachindex(eqs)
        eqs[i] = substitute_ẋ_with_v(eqs[i], ẋ, v)
    end
end

function substitute_lagrangian_variables(equ, x, v)
    @variables X[axes(x, 1)]
    @variables V[axes(v, 1)]
    substitute(equ, [z=>Z for (z,Z) in zip([x..., v...], [X..., V...])])
end

function substitute_lagrangian_variables!(eqs, x, v)
    for i in eachindex(eqs)
        eqs[i] = substitute_lagrangian_variables(eqs[i], x, v)
    end
end


struct LagrangianSystem
    L
    t
    x
    v
    parameters
    equations
    functions

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
            substitute_lagrangian_variables!(eq, x, v)
        end

        L = substitute_lagrangian_variables(L, x, v)
        a = inv(M) * (f - N * V)
        ϕ = P .- ϑ
        ψ = F .- g

        equs = (
            L = L,
            EL = EL,
            a = a,
            f = f,
            g = g,
            ϑ = ϑ,
            θ = θ,
            ω = ω,
            Ω = Ω,
            ϕ = ϕ,
            ψ = ψ,
            M = M,
            N = N,
            Σ = Σ,
        )

        code = (
            L  = substitute_parameters(build_function(equs.L,  t, X, V, params...), params),
            EL = substitute_parameters(build_function(equs.EL, t, X, V, params...)[2], params),
            a  = substitute_parameters(build_function(equs.a,  t, X, V, params...)[2], params),
            f  = substitute_parameters(build_function(equs.f,  t, X, V, params...)[2], params),
            g  = substitute_parameters(build_function(equs.g,  t, X, V, params...)[2], params),
            p  = substitute_parameters(build_function(equs.ϑ,  t, X, V, params...)[1], params),
            ϑ  = substitute_parameters(build_function(equs.ϑ,  t, X, V, params...)[2], params),
            θ  = substitute_parameters(build_function(equs.θ,  t, X, V, params...)[2], params),
            ω  = substitute_parameters(build_function(equs.ω,  t, X, V, params...)[2], params),
            Ω  = substitute_parameters(build_function(equs.Ω,  t, X, V, params...)[2], params),
            ϕ  = substitute_parameters(build_function(equs.ϕ,  t, X, V, P, params...)[2], params),
            ψ  = substitute_parameters(build_function(equs.ψ,  t, X, V, P, F, params...)[2], params),
            M  = substitute_parameters(build_function(equs.M,  t, X, V, params...)[2], params),
            P  = substitute_parameters(build_function(equs.Σ,  t, X, V, params...)[2], params),
        )

        funcs = NamedTuple{keys(code)}(Tuple(@RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(c)) for c in code))
        # (
        #     EL = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.EL)),
        #     a  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.a)),
        #     f  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.f)),
        #     g  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.g)),
        #     p  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.p)),
        #     ϑ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.ϑ)),
        #     ω  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.ω)),
        #     Ω  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.Ω)),
        #     ϕ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.ϕ)),
        #     ψ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.ψ)),
        #     L  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.L)),
        #     M  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.M)),
        #     P  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code.P)),
        # )

        funcs_param = length(params) > 0 ? funcs : (
            L = (t,x,v,params)     -> funcs.L(t,x,v),
            EL= (e,t,x,v,params)   -> funcs.EL(e,t,x,v),
            a = (a,t,x,v,params)   -> funcs.a(a,t,x,v),
            ϑ = (ϑ,t,x,v,params)   -> funcs.ϑ(ϑ,t,x,v),
            θ = (θ,t,x,v,params)   -> funcs.θ(θ,t,x,v),
            f = (f,t,x,v,params)   -> funcs.f(f,t,x,v),
            g = (f̄,t,x,v,λ,params) -> funcs.f̄(f̄,t,x,λ),
            ω = (ω,t,x,v,params)   -> funcs.ω(ω,t,x,v),
            Ω = (Ω,t,x,v,params)   -> funcs.Ω(Ω,t,x,v),
            ϕ = (ϕ,t,x,v,params)   -> funcs.ϕ(ϕ,t,x,v),
            ψ = (ψ,t,x,v,params)   -> funcs.ψ(ψ,t,x,v),
            v̄ = (v,t,x,p,params)   -> funcs.ẋ(v,t,x,p),
            f̄ = (f,t,x,v,params)   -> funcs.f(f,t,x,v),
            M = (M,t,x,v,params)   -> funcs.M(M,t,x,v),
            P = (P,t,x,v,params)   -> funcs.P(P,t,x,v),
        )

        new(L, t, x, v, params, equs, funcs_param)
    end
end

lagrangian(lsys::LagrangianSystem) = lsys.L
parameters(lsys::LagrangianSystem) = lsys.parameters
variables(lsys::LagrangianSystem) = (lsys.t, lsys.x, lsys.v)
equations(lsys::LagrangianSystem) = lsys.equations
functions(lsys::LagrangianSystem) = lsys.functions


function lagrangian_variables(dimension::Int)
    t = parameter(:t)
    
    @variables (x(t))[1:dimension]
    @variables (v(t))[1:dimension]

    return (t,x,v)
end


function LODE(lsys::LagrangianSystem)
    eqs = functions(lsys)
    LODE(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L; v̄ = eqs.v̄, f̄ = eqs.f̄)
end

function LODEProblem(lsys::LagrangianSystem, tspan::Tuple, tstep::Real, ics::NamedTuple)
    eqs = functions(lsys)
    LODEProblem(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, tspan, tstep, ics; v̄ = eqs.v̄, f̄ = eqs.f̄)
end

function LODEProblem(lsys::LagrangianSystem, tspan::Tuple, tstep::Real, q₀::State, p₀::State, λ₀::State = zero(q₀))
    ics = (q = q₀, p = p₀, λ = λ₀)
    LODEProblem(lsys, tspan, tstep, ics)
end
