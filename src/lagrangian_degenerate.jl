"""
    DegenerateLagrangianSystem
"""
struct DegenerateLagrangianSystem
    L
    t
    x
    v
    parameters
    equations
    functions

    function DegenerateLagrangianSystem(θ, H, t, x, v, params = NamedTuple())

        @assert eachindex(x) == eachindex(v)

        @variables (p(t))[axes(x, 1)]
        @variables (f(t))[axes(x, 1)]
    
        @variables X[axes(x, 1)]
        @variables V[axes(v, 1)]
        @variables P[axes(v, 1)]
        @variables F[axes(x, 1)]
        @variables Λ[axes(v, 1)]

        Dt = Differential(t)
        Dx = collect(Differential.(x))
        Dv = collect(Differential.(v))
        ẋ  = collect(Dt.(x))

        L  = θ ⋅ v - H
        EL = [expand_derivatives(Dx[i](L) - Dt(Dv[i](L))) for i in eachindex(Dx,Dv)]
        ∇H = [expand_derivatives(dx(H)) for dx in Dx]
        ϑ  = [expand_derivatives(dv(L)) for dv in Dv]
        f  = [expand_derivatives(dx(L)) for dx in Dx]
        g  = [expand_derivatives(Dt.(ϑ))]
        ω  = [expand_derivatives(simplify(Dx[i](ϑ[j]) - Dx[j](ϑ[i]))) for i in eachindex(Dx,ϑ), j in eachindex(Dx,ϑ)]
        N  = [expand_derivatives(simplify(Dx[i](ϑ[j]))) for i in eachindex(Dv), j in eachindex(ϑ)]

        equs = (
            L = L,
            EL = EL,
            ∇H = ∇H,
            f = f,
            g = g,
            ϑ = ϑ,
            ω = ω,
            N = N,
        )

        equs_subs = substitute_lagrangian_variables(equs, x, ẋ, v)
        equs_subs = merge(equs_subs, (
            ϕ = P .- equs_subs.ϑ,
            ψ = F .- equs_subs.g,
            σ = simplify.(inv(equs_subs.ω)),
        ))

        equs_subs = merge(equs_subs, (
            ẋ = equs_subs.σ * equs_subs.∇H,
        ))

        # equs = substitute_v_with_ẋ(equs, v, ẋ)
        # equs = merge(equs, (
        #     ϕ = p .- equs.ϑ,
        #     ψ = ṗ .- equs.g,
        # ))

        code = (
            L  = substitute_parameters(build_function(equs_subs.L,  t, X, V, params...; nanmath = false), params),
            EL = substitute_parameters(build_function(equs_subs.EL, t, X, V, params...; nanmath = false)[2], params),
            ∇H = substitute_parameters(build_function(equs_subs.∇H, t, X, V, params...; nanmath = false)[2], params),
            ẋ  = substitute_parameters(build_function(equs_subs.ẋ,  t, X, V, params...; nanmath = false)[2], params),
            f  = substitute_parameters(build_function(equs_subs.f,  t, X, V, params...; nanmath = false)[2], params),
            g  = substitute_parameters(build_function(equs_subs.g,  t, X, Λ, V, params...; nanmath = false)[2], params),
            p  = substitute_parameters(build_function(equs_subs.ϑ,  t, X, V, params...; nanmath = false)[1], params),
            ϑ  = substitute_parameters(build_function(equs_subs.ϑ,  t, X, V, params...; nanmath = false)[2], params),
            ω  = substitute_parameters(build_function(equs_subs.ω,  t, X, V, params...; nanmath = false)[2], params),
            ϕ  = substitute_parameters(build_function(equs_subs.ϕ,  t, X, V, P, params...; nanmath = false)[2], params),
            ψ  = substitute_parameters(build_function(equs_subs.ψ,  t, X, V, P, F, params...; nanmath = false)[2], params),
            P  = substitute_parameters(build_function(equs_subs.σ,  t, X, V, params...; nanmath = false)[2], params),
        )

        funcs = generate_code(code)

        new(L, t, x, v, params, equs, funcs)
    end
end

lagrangian(lsys::DegenerateLagrangianSystem) = lsys.L
parameters(lsys::DegenerateLagrangianSystem) = lsys.parameters
variables(lsys::DegenerateLagrangianSystem) = (lsys.t, lsys.x, lsys.v)
equations(lsys::DegenerateLagrangianSystem) = lsys.equations
functions(lsys::DegenerateLagrangianSystem) = lsys.functions

function Base.show(io::IO, lsys::DegenerateLagrangianSystem)
    print(io, "\nDegenerate Lagrangian system with\n")
    print(io, "\nL = ")
    print(io, lagrangian(lsys))
    print(io, "\n\nand equations of motion\n\n")
    for eq in equations(lsys).EL
        print(io, eq)
        print(io, "\n")
    end
end


function LODE(lsys::DegenerateLagrangianSystem; v̄ = functions(lsys).ẋ, f̄ = functions(lsys).f, kwargs...)
    eqs = functions(lsys)
    LODE(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L; v̄ = v̄, f̄ = f̄, kwargs...)
end

function LODEProblem(lsys::DegenerateLagrangianSystem, tspan::Tuple, tstep::Real, ics::NamedTuple; v̄ = functions(lsys).ẋ, f̄ = functions(lsys).f, kwargs...)
    eqs = functions(lsys)
    LODEProblem(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, tspan, tstep, ics; v̄ = v̄, f̄ = f̄, kwargs...)
end

function LODEProblem(lsys::DegenerateLagrangianSystem, tspan::Tuple, tstep::Real, q₀::StateVariable, p₀::StateVariable, λ₀::AlgebraicVariable; kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    LODEProblem(lsys, tspan, tstep, ics; kwargs...)
end

function LODEProblem(lsys::DegenerateLagrangianSystem, tspan::Tuple, tstep::Real, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray = zero(q₀); kwargs...)
    _q₀ = StateVariable(q₀)
    _p₀ = StateVariable(p₀)
    _λ₀ = AlgebraicVariable(λ₀)
    LODEProblem(lsys, tspan, tstep, _q₀, _p₀, _λ₀; kwargs...)
end
