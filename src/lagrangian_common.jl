
function substitute_v_with_ẋ(equ, v, ẋ)
    for i in eachindex(ẋ, v)
        equ = substitute(equ, v[i] => ẋ[i])
    end
    return equ
end

function substitute_v_with_ẋ(equs::AbstractArray, v, ẋ)
    [substitute_v_with_ẋ(eq, v, ẋ) for eq in equs]
end

function substitute_v_with_ẋ(equs::NamedTuple, v, ẋ)
    NamedTuple{keys(equs)}(Tuple(substitute_v_with_ẋ(eq, v, ẋ) for eq in equs))
end

function substitute_ẋ_with_v(equ, ẋ, v)
    substitute(equ, [ẋᵢ => vᵢ for (ẋᵢ, vᵢ) in zip(ẋ, v)])
end

function substitute_ẋ_with_v(equs::AbstractArray, ẋ, v)
    [substitute_ẋ_with_v(eq, ẋ, v) for eq in equs]
end

"""
    substitute_acceleration(equ, dv, Λ)

Replace the time derivatives of the velocities, `dv = d/dt(v)`, by the algebraic variable `Λ`.

Any quantity built from `Dt(∂L/∂v)` — `g`, and hence `EL = f - g` and `ψ = ṗ - g` — contains the
acceleration. `substitute_ẋ_with_v` removes only the *first* derivatives `d/dt(x)`, so without this
the generated code refers to an unevaluated `Differential(t)(v[i])`: a variable that does not exist,
which makes the function throw as soon as it is called. `Λ` is the slot the generated signatures
already carry for exactly this purpose.
"""
function substitute_acceleration(equ, dv, Λ)
    substitute(equ, [dvᵢ => Λᵢ for (dvᵢ, Λᵢ) in zip(dv, collect(Λ))])
end

function substitute_acceleration(equs::AbstractArray, dv, Λ)
    [substitute_acceleration(eq, dv, Λ) for eq in equs]
end

function substitute_acceleration(equs::NamedTuple, dv, Λ)
    NamedTuple{keys(equs)}(Tuple(substitute_acceleration(eq, dv, Λ) for eq in equs))
end

function substitute_lagrangian_variables(equ, x, v)
    @variables X[axes(x, 1)]
    @variables V[axes(v, 1)]
    substitute(equ, [zᵢ => Zᵢ for (zᵢ, Zᵢ) in zip([x..., v...], [X..., V...])])
end

function substitute_lagrangian_variables(equs::AbstractArray, x, v)
    [substitute_lagrangian_variables(eq, x, v) for eq in equs]
end

function substitute_lagrangian_variables(equs::NamedTuple, x, ẋ, v)
    NamedTuple{keys(equs)}(Tuple(substitute_lagrangian_variables(substitute_ẋ_with_v(eq, ẋ, v), x, v)
    for eq in equs))
end

function lagrangian_variables(dimension::Int)
    @variables t
    @variables x(t)[1:dimension]
    @variables v(t)[1:dimension]

    return (t, x, v)
end

function lagrangian_derivatives(t, x, v)
    Dt = Differential(t)
    Dx = collect(Differential.(x))
    Dv = collect(Differential.(v))

    return (Dt, Dx, Dv)
end
