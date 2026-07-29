
"""
    substitute_parameters(code, params)

Rewrite the anonymous function `code` built by `Symbolics.build_function` so that it takes a
single `params` argument in place of one argument per entry of `params`.

`symbolize` names the symbolic stand-in for parameter `k` as `kₚ`, so `build_function` emits a
signature like `(ˍ₋out, t, X, V, Gₚ, mₚ)` and a body referring to `Gₚ` and `mₚ`. This returns
the equivalent `(ˍ₋out, t, X, V, params)` with `params.G` and `params.m` in the body. When
`params` is empty, only the `params` argument is appended.

This walks the expression rather than round-tripping it through `string` and `Meta.parse`. The
generated code reaches ~1 MB for a matrix-valued quantity in 18 degrees of freedom, where the
text round-trip is both slow and dependent on how `Base.show` happens to render an `Expr`.
"""
function substitute_parameters(code, params)
    Meta.isexpr(code, :function, 2) && Meta.isexpr(code.args[1], :tuple) ||
        error("substitute_parameters expects `function (args...) ... end`, got $(repr(code))")

    substitutions = Dict{Symbol,Expr}(Symbol("$(k)ₚ") => :(params.$k) for k in keys(params))

    arguments = filter(a -> !(a isa Symbol && haskey(substitutions, a)), code.args[1].args)
    signature = Expr(:tuple, arguments..., :params)

    return Expr(:function, signature, substitute_symbols(code.args[2], substitutions))
end

substitute_symbols(ex::Symbol, substitutions) = get(substitutions, ex, ex)

function substitute_symbols(ex::Expr, substitutions)
    Expr(ex.head, Any[substitute_symbols(a, substitutions) for a in ex.args]...)
end

# Literals, `LineNumberNode`s and `QuoteNode`s pass through unchanged.
substitute_symbols(ex, substitutions) = ex


function symbolize(p::Union{AbstractArray,Tuple}, name)
    vars = @variables $(name)[axes(p)...]
    first(vars)
end

function symbolize(x::T, name) where {T<:Number}
    vars = @variables $(name)::Real
    first(vars)
end

function symbolize(p::Union{T,AbstractArray{T},Symbolics.Arr{<:Num}}, name) where {T<:SymbolicUtils.BasicSymbolicImpl.Type}
    p
end

function symbolize(params::NamedTuple)
    NamedTuple{keys(params)}(Tuple(symbolize(v, Symbol("$(k)ₚ")) for (k, v) in pairs(params)))
end

function symbolize(::Union{Nothing,NullParameters})
    NamedTuple{}()
end


function generate_code(code)
    @RuntimeGeneratedFunction(code)
end

function generate_code(code::NamedTuple)
    NamedTuple{keys(code)}(Tuple(generate_code(c) for c in code))
end
