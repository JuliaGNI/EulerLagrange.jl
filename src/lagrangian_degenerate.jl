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

    function DegenerateLagrangianSystem(K, H, t, x, v, params=NamedTuple(); simplify=true, scalarize=true)

        @assert eachindex(x) == eachindex(v)

        @variables p(t)[axes(x, 1)]
        @variables f(t)[axes(x, 1)]

        @variables X[axes(x, 1)]
        @variables V[axes(v, 1)]
        @variables P[axes(v, 1)]
        @variables F[axes(x, 1)]
        @variables Λ[axes(v, 1)]

        Dt, Dx, Dv = lagrangian_derivatives(t, x, v)

        ẋ = collect(Dt.(x))

        Ks = scalarize ? Symbolics.scalarize(K) : K
        Hs = scalarize ? Symbolics.scalarize(H) : H

        Ks = simplify ? Symbolics.simplify(Ks) : Ks
        Hs = simplify ? Symbolics.simplify(Hs) : Hs

        Ls = Ks - Hs
        ∇H = [expand_derivatives(dx(Hs)) for dx in Dx]
        ϑ = [expand_derivatives(dv(Ls)) for dv in Dv]
        f = [expand_derivatives(dx(Ls)) for dx in Dx]
        g = [expand_derivatives(dx(Ks)) for dx in Dx]
        ḡ = [expand_derivatives(Dt(θ)) for θ in ϑ]
        ω = [expand_derivatives(Symbolics.simplify(Dx[i](ϑ[j]) - Dx[j](ϑ[i]))) for i in eachindex(Dx, ϑ), j in eachindex(Dx, ϑ)]
        N = [expand_derivatives(Symbolics.simplify(Dx[i](ϑ[j]))) for i in eachindex(Dv), j in eachindex(ϑ)]
        u = [u for u in ẋ]
        ū = [u for u in ẋ]
        EL = [f[i] - g[i] for i in eachindex(f, g)]

        equs = (
            L=Ls,
            H=Hs,
            EL=EL,
            ∇H=∇H,
            f=f,
            u=u,
            g=g,
            ū=ū,
            ḡ=ḡ,
            ϑ=ϑ,
            ω=ω,
            N=N,
        )

        equs_subs = substitute_lagrangian_variables(equs, x, ẋ, v)

        σ = inv(equs_subs.ω)

        equs_subs = merge(equs_subs, (
            ϕ=[P[i] - equs_subs.ϑ[i] for i in eachindex(P, equs_subs.ϑ)],
            ψ=[F[i] - equs_subs.ḡ[i] for i in eachindex(F, equs_subs.ḡ)],
            σ=simplify ? Symbolics.simplify.(σ) : σ,
        ))

        ẋeq = equs_subs.σ * equs_subs.∇H

        equs_subs = merge(equs_subs, (
            ẋ=simplify ? Symbolics.simplify.(ẋeq) : ẋeq,
        ))

        code = (
            L=substitute_parameters(build_function(equs_subs.L, t, X, V, params...; nanmath=false), params),
            H=substitute_parameters(build_function(equs_subs.H, t, X, params...; nanmath=false), params),
            EL=substitute_parameters(build_function(equs_subs.EL, t, X, V, params...; nanmath=false)[2], params),
            ∇H=substitute_parameters(build_function(equs_subs.∇H, t, X, V, params...; nanmath=false)[2], params),
            ẋ=substitute_parameters(build_function(equs_subs.ẋ, t, X, params...; nanmath=false)[2], params),
            v=substitute_parameters(build_function(equs_subs.ẋ, t, X, P, params...; nanmath=false)[2], params),
            f=substitute_parameters(build_function(equs_subs.f, t, X, V, params...; nanmath=false)[2], params),
            u=substitute_parameters(build_function(equs_subs.u, t, X, Λ, V, params...; nanmath=false)[2], params),
            g=substitute_parameters(build_function(equs_subs.g, t, X, Λ, V, params...; nanmath=false)[2], params),
            ū=substitute_parameters(build_function(equs_subs.ū, t, X, Λ, P, V, params...; nanmath=false)[2], params),
            ḡ=substitute_parameters(build_function(equs_subs.ḡ, t, X, Λ, P, V, params...; nanmath=false)[2], params),
            p=substitute_parameters(build_function(equs_subs.ϑ, t, X, V, params...; nanmath=false)[1], params),
            ϑ=substitute_parameters(build_function(equs_subs.ϑ, t, X, V, params...; nanmath=false)[2], params),
            ω=substitute_parameters(build_function(equs_subs.ω, t, X, V, params...; nanmath=false)[2], params),
            ϕ=substitute_parameters(build_function(equs_subs.ϕ, t, X, V, P, params...; nanmath=false)[2], params),
            ψ=substitute_parameters(build_function(equs_subs.ψ, t, X, V, P, F, params...; nanmath=false)[2], params),
            P=substitute_parameters(build_function(equs_subs.σ, t, X, V, params...; nanmath=false)[2], params),
        )

        new(Ls, t, x, v, params, equs, generate_code(code))
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
    # print(io, "\n\nand equations of motion\n\n")
    # for eq in equations(lsys).EL
    #     print(io, eq)
    #     print(io, "\n")
    # end
end


function ODE(lsys::DegenerateLagrangianSystem; kwargs...)
    eqs = functions(lsys)
    ODE(eqs.ẋ; invariants=(h=eqs.H,), kwargs...)
end

function ODEProblem(lsys::DegenerateLagrangianSystem, tspan::Tuple, tstep::Real, ics...; kwargs...)
    eqs = functions(lsys)
    ODEProblem(eqs.ẋ, tspan, tstep, ics...; invariants=(h=eqs.H,), kwargs...)
end


function LODE(lsys::DegenerateLagrangianSystem; v̄=functions(lsys).v, f̄=functions(lsys).f, kwargs...)
    eqs = functions(lsys)
    LODE(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L; v̄=v̄, f̄=f̄, invariants=(h=(t, q, v, params) -> eqs.H(t, q, params),), kwargs...)
end

function LODEProblem(lsys::DegenerateLagrangianSystem, tspan::Tuple, tstep::Real, ics...; v̄=functions(lsys).v, f̄=functions(lsys).f, kwargs...)
    eqs = functions(lsys)
    LODEProblem(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, tspan, tstep, ics...; v̄=v̄, f̄=f̄, invariants=(h=(t, q, v, params) -> eqs.H(t, q, params),), kwargs...)
end


function LDAE(lsys::DegenerateLagrangianSystem; v̄=functions(lsys).v, f̄=functions(lsys).f, kwargs...)
    eqs = functions(lsys)
    LDAE(eqs.ϑ, eqs.f, (u, t, q, v, p, λ, params) -> eqs.u(u, t, q, v, λ, params), (g, t, q, v, p, λ, params) -> eqs.g(g, t, q, v, λ, params), eqs.ϕ, eqs.ū, eqs.ḡ, eqs.ψ, eqs.ω, eqs.L; v̄=v̄, f̄=f̄, invariants=(h=(t, q, v, params) -> eqs.H(t, q, params),), kwargs...)
end

function LDAEProblem(lsys::DegenerateLagrangianSystem, tspan::Tuple, tstep::Real, ics...; v̄=functions(lsys).v, f̄=functions(lsys).f, kwargs...)
    eqs = functions(lsys)
    LDAEProblem(eqs.ϑ, eqs.f, (u, t, q, v, p, λ, params) -> eqs.u(u, t, q, v, λ, params), (g, t, q, v, p, λ, params) -> eqs.g(g, t, q, v, λ, params), eqs.ϕ, eqs.ū, eqs.ḡ, eqs.ψ, eqs.ω, eqs.L, tspan, tstep, ics...; v̄=v̄, f̄=f̄, invariants=(h=(t, q, v, params) -> eqs.H(t, q, params),), kwargs...)
end
