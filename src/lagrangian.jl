
function substitute_v_with_ẋ(equ, v, ẋ)
    for i in eachindex(ẋ,v)
        equ = substitute(equ, v[i] => ẋ[i])
    end
    return equ
end

function substitute_v_with_ẋ(equs::Union{AbstractArray, ArrayLike}, v, ẋ)
    [substitute_v_with_ẋ(eq, v, ẋ) for eq in equs]
end

function substitute_v_with_ẋ(equs::NamedTuple, v, ẋ)
    NamedTuple{keys(equs)}(Tuple(substitute_v_with_ẋ(eq, v, ẋ) for eq in equs))
end

function substitute_ẋ_with_v(equ, ẋ, v)
    substitute(equ, [ẋᵢ=>vᵢ for (ẋᵢ,vᵢ) in zip(ẋ,v)])
end

function substitute_ẋ_with_v(equs::Union{AbstractArray, ArrayLike}, ẋ, v)
    [substitute_ẋ_with_v(eq, ẋ, v) for eq in equs]
end

function substitute_lagrangian_variables(equ, x, v)
    @variables X[axes(x, 1)]
    @variables V[axes(v, 1)]
    substitute(equ, [zᵢ=>Zᵢ for (zᵢ,Zᵢ) in zip([x..., v...], [X..., V...])])
end

function substitute_lagrangian_variables(equs::Union{AbstractArray, ArrayLike}, x, v)
    [substitute_lagrangian_variables(eq, x, v) for eq in equs]
end

function substitute_lagrangian_variables(equs::NamedTuple, x, ẋ, v)
    NamedTuple{keys(equs)}(Tuple(substitute_lagrangian_variables(substitute_ẋ_with_v(eq, ẋ, v), x, v) for eq in equs))
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

        @variables (p(t))[axes(x, 1)]
        @variables (f(t))[axes(x, 1)]
    
        @variables X[axes(x, 1)]
        @variables V[axes(v, 1)]
        @variables P[axes(v, 1)]
        @variables F[axes(x, 1)]

        Dt = Differential(t)
        Dx = collect(Differential.(x))
        Dv = collect(Differential.(v))
        Dz = vcat(Dx,Dv)
        ẋ  = collect(Dt.(x))
        ṗ  = collect(Dt.(p))

        EL = [expand_derivatives(Dx[i](L) - Dt(Dv[i](L))) for i in eachindex(Dx,Dv)]
        f  = [expand_derivatives(dx(L)) for dx in Dx]
        g  = [expand_derivatives(Dt(dv(L))) for dv in Dv]
        ϑ  = [expand_derivatives(dv(L)) for dv in Dv]
        θ  = [expand_derivatives(dz(L)) for dz in Dz]
        Ω  = [expand_derivatives(simplify(Dx[i](ϑ[j]) - Dx[j](ϑ[i]))) for i in eachindex(Dx,ϑ), j in eachindex(Dx,ϑ)]
        ω  = [expand_derivatives(simplify(Dz[i](θ[j]) - Dz[j](θ[i]))) for i in eachindex(Dz,θ), j in eachindex(Dz,θ)]
        M  = [expand_derivatives(simplify(Dv[i](ϑ[j]))) for i in eachindex(Dv), j in eachindex(ϑ)]
        N  = [expand_derivatives(simplify(Dx[i](ϑ[j]))) for i in eachindex(Dv), j in eachindex(ϑ)]

        equs = (
            L = L,
            EL = EL,
            f = f,
            g = g,
            ϑ = ϑ,
            θ = θ,
            Ω = Ω,
            ω = ω,
            M = M,
            N = N,
        )

        equs_subs = substitute_lagrangian_variables(equs, x, ẋ, v)
        equs_subs = merge(equs_subs, (
            a = inv(equs_subs.M) * (equs_subs.f - equs_subs.N * V),
            ϕ = P .- equs_subs.ϑ,
            ψ = F .- equs_subs.g,
            σ = simplify.(inv(equs_subs.ω)),
            Σ = simplify.(inv(equs_subs.Ω)),
        ))

        equs = substitute_v_with_ẋ(equs, v, ẋ)
        equs = merge(equs, (
            ϕ = p .- equs.ϑ,
            ψ = ṗ .- equs.g,
        ))

        code = (
            L  = substitute_parameters(build_function(equs_subs.L,  t, X, V, params...), params),
            EL = substitute_parameters(build_function(equs_subs.EL, t, X, V, params...)[2], params),
            a  = substitute_parameters(build_function(equs_subs.a,  t, X, V, params...)[2], params),
            f  = substitute_parameters(build_function(equs_subs.f,  t, X, V, params...)[2], params),
            g  = substitute_parameters(build_function(equs_subs.g,  t, X, V, params...)[2], params),
            p  = substitute_parameters(build_function(equs_subs.ϑ,  t, X, V, params...)[1], params),
            ϑ  = substitute_parameters(build_function(equs_subs.ϑ,  t, X, V, params...)[2], params),
            θ  = substitute_parameters(build_function(equs_subs.θ,  t, X, V, params...)[2], params),
            ω  = substitute_parameters(build_function(equs_subs.ω,  t, X, V, params...)[2], params),
            Ω  = substitute_parameters(build_function(equs_subs.Ω,  t, X, V, params...)[2], params),
            ϕ  = substitute_parameters(build_function(equs_subs.ϕ,  t, X, V, P, params...)[2], params),
            ψ  = substitute_parameters(build_function(equs_subs.ψ,  t, X, V, P, F, params...)[2], params),
            M  = substitute_parameters(build_function(equs_subs.M,  t, X, V, params...)[2], params),
            P  = substitute_parameters(build_function(equs_subs.Σ,  t, X, V, params...)[2], params),
        )

        funcs = generate_code(code)

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

function Base.show(io::IO, lsys::LagrangianSystem)
    print(io, "\nLagrangian system with\n")
    print(io, "\nL = ")
    print(io, lagrangian(lsys))
    print(io, "\n\nand equations of motion\n\n")
    for eq in equations(lsys).EL
        print(io, eq)
        print(io, "\n")
    end
end

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