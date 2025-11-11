
function substitute_v_with_ẋ(equ, v, ẋ)
    for i in eachindex(ẋ, v)
        equ = substitute(equ, v[i] => ẋ[i])
    end
    return equ
end

function substitute_v_with_ẋ(equs::AbstractArray, v, ẋ)
    [substitute_v_with_ẋ(eq, v, ẋ) for eq in equs]
end

function substitute_v_with_ẋ(equs::NamedTuple, v, ẋ)
    NamedTuple{keys(equs)}(Tuple(substitute_v_with_ẋ(eq, v, ẋ) for eq in equs))
end

function substitute_ẋ_with_v(equ, ẋ, v)
    substitute(equ, [ẋᵢ => vᵢ for (ẋᵢ, vᵢ) in zip(ẋ, v)])
end

function substitute_ẋ_with_v(equs::AbstractArray, ẋ, v)
    [substitute_ẋ_with_v(eq, ẋ, v) for eq in equs]
end

function substitute_lagrangian_variables(equ, x, v)
    @variables X[axes(x, 1)]
    @variables V[axes(v, 1)]
    substitute(equ, [zᵢ => Zᵢ for (zᵢ, Zᵢ) in zip([x..., v...], [X..., V...])])
end

function substitute_lagrangian_variables(equs::AbstractArray, x, v)
    [substitute_lagrangian_variables(eq, x, v) for eq in equs]
end

function substitute_lagrangian_variables(equs::NamedTuple, x, ẋ, v)
    NamedTuple{keys(equs)}(Tuple(substitute_lagrangian_variables(substitute_ẋ_with_v(eq, ẋ, v), x, v) for eq in equs))
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
