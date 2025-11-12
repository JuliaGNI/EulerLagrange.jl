"""
    LagrangianSystem
"""
struct LagrangianSystem
    L
    t
    x
    v
    parameters
    equations
    functions

    function LagrangianSystem(L, t, x, v, params=NamedTuple(); simplify=true, scalarize=true)

        @assert eachindex(x) == eachindex(v)

        @variables p(t)[axes(x, 1)]
        @variables f(t)[axes(x, 1)]

        @variables X[axes(x, 1)]
        @variables V[axes(v, 1)]
        @variables P[axes(v, 1)]
        @variables F[axes(x, 1)]
        @variables Λ[axes(v, 1)]

        Dt, Dx, Dv = lagrangian_derivatives(t, x, v)

        Dz = vcat(Dx, Dv)
        ẋ = collect(Dt.(x))
        ṗ = collect(Dt.(p))

        Ls = scalarize ? Symbolics.scalarize(L) : L
        Ls = simplify ? Symbolics.simplify(Ls) : Ls
        Ls = Num(Ls)

        f = [expand_derivatives(dx(Ls)) for dx in Dx]
        g = [expand_derivatives(Dt(dv(Ls))) for dv in Dv]
        ϑ = [expand_derivatives(dv(Ls)) for dv in Dv]
        θ = [expand_derivatives(dz(Ls)) for dz in Dz]
        EL = [f[i] - g[i] for i in eachindex(f, g)]

        Ω = [Dx[i](ϑ[j]) - Dx[j](ϑ[i]) for i in eachindex(Dx, ϑ), j in eachindex(Dx, ϑ)]
        ω = [Dz[i](θ[j]) - Dz[j](θ[i]) for i in eachindex(Dz, θ), j in eachindex(Dz, θ)]
        M = [Dv[i](ϑ[j]) for i in eachindex(Dv), j in eachindex(ϑ)]
        N = [Dx[i](ϑ[j]) for i in eachindex(Dv), j in eachindex(ϑ)]

        if simplify
            Ω = Symbolics.simplify.(Ω)
            ω = Symbolics.simplify.(ω)
            M = Symbolics.simplify.(M)
            N = Symbolics.simplify.(N)
        end

        Ω = expand_derivatives.(Ω)
        ω = expand_derivatives.(ω)
        M = expand_derivatives.(M)
        N = expand_derivatives.(N)

        equs = (
            L=Ls,
            EL=EL,
            f=f,
            g=g,
            ϑ=ϑ,
            θ=θ,
            Ω=Ω,
            ω=ω,
            M=M,
            N=N,
        )

        # _simplify(expr, dosimplify) = dosimplify ? simplify.(expr) : expr

        equs_subs = substitute_lagrangian_variables(equs, x, ẋ, v)
        equs_subs = merge(equs_subs, (
            # a = inv(equs_subs.M) * (equs_subs.f - equs_subs.N * V),
            ϕ=P .- equs_subs.ϑ,
            ψ=F .- equs_subs.g,
            # σ = _simplify(inv(equs_subs.ω), dosimplify),
            # Σ = _simplify(inv(equs_subs.Ω), dosimplify),
        ))

        equs = substitute_v_with_ẋ(equs, v, ẋ)
        equs = merge(equs, (
            ϕ=p .- equs.ϑ,
            ψ=ṗ .- equs.g,
        ))

        code = (
            L=substitute_parameters(build_function(equs_subs.L, t, X, V, params...; nanmath=false), params),
            EL=substitute_parameters(build_function(equs_subs.EL, t, X, V, params...; nanmath=false)[2], params),
            # a  = substitute_parameters(build_function(equs_subs.a,  t, X, V, params...; nanmath = false)[2], params),
            f=substitute_parameters(build_function(equs_subs.f, t, X, V, params...; nanmath=false)[2], params),
            g=substitute_parameters(build_function(equs_subs.g, t, X, Λ, V, params...; nanmath=false)[2], params),
            p=substitute_parameters(build_function(equs_subs.ϑ, t, X, V, params...; nanmath=false)[1], params),
            ϑ=substitute_parameters(build_function(equs_subs.ϑ, t, X, V, params...; nanmath=false)[2], params),
            θ=substitute_parameters(build_function(equs_subs.θ, t, X, V, params...; nanmath=false)[2], params),
            ω=substitute_parameters(build_function(equs_subs.ω, t, X, V, params...; nanmath=false)[2], params),
            Ω=substitute_parameters(build_function(equs_subs.Ω, t, X, V, params...; nanmath=false)[2], params),
            ϕ=substitute_parameters(build_function(equs_subs.ϕ, t, X, V, P, params...; nanmath=false)[2], params),
            ψ=substitute_parameters(build_function(equs_subs.ψ, t, X, V, P, F, params...; nanmath=false)[2], params),
            M=substitute_parameters(build_function(equs_subs.M, t, X, V, params...; nanmath=false)[2], params),
            # P  = substitute_parameters(build_function(equs_subs.Σ,  t, X, V, params...; nanmath = false)[2], params),
        )

        new(Ls, t, x, v, params, equs, generate_code(code))
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
    # print(io, "\n\nand equations of motion\n\n")
    # for eq in equations(lsys).EL
    #     print(io, eq)
    #     print(io, "\n")
    # end
end


function LODE(lsys::LagrangianSystem; kwargs...)
    eqs = functions(lsys)
    LODE(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L; kwargs...)
end

function LODEProblem(lsys::LagrangianSystem, tspan::Tuple, tstep::Real, ics...; kwargs...)
    eqs = functions(lsys)
    LODEProblem(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, tspan, tstep, ics...; kwargs...)
end
