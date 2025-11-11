
function substitute_parameters(code, params)
    if length(params) > 0
        # generate string of parameter arguments as they appear in generated code
        paramstr = string(Tuple(string(k) * "ₚ, " for k in keys(params))...)
        # remove trailing comma
        paramstr = paramstr[begin:prevind(paramstr, first((findlast(",", paramstr))))]
        # convert code expression to string
        code_str = string(code)
        # extract function header
        func_str = code_str[first(findfirst("function (", code_str)):last(findfirst(paramstr * ")", code_str))]
        # replace parameter list in function header with `params`
        func_str_params = replace(func_str, paramstr * ")" => "params)")
        # replace function header in code
        code_str = replace(code_str, func_str => func_str_params)
        # replace all params with named tuple entries
        for k in keys(params)
            code_str = replace(code_str, "$(k)ₚ" => "params.$(k)")
        end
        # convert code string back to expression
        return Meta.parse(code_str)
    else
        # convert code expression to string
        code_str = string(code)
        # extract function header
        func_str = code_str[first(findfirst("function (", code_str)):last(findfirst(")", code_str))]
        # append params argument to function header
        func_str_params = replace(func_str, ")" => ", params)")
        # replace function header in code
        code_str = replace(code_str, func_str => func_str_params)
        # convert code string back to expression
        return Meta.parse(code_str)
    end
end


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
