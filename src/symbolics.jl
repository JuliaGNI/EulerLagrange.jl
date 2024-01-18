
parameter(name::Symbol) = Num(Sym{Real}(name))


function substitute_parameters(code, params)
    if length(params) > 0
        # generate string of parameter arguments as they appear in generated code
        paramstr = string(Tuple("params_" * string(k) * ", " for k in keys(params))...)
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
        # replace all remaining params with named tuple entries
        code_str = replace(code_str, "params_" => "params.")
        # convert code string back to expression
        return Meta.parse(code_str)
    else
        return code
    end
end


function symbolize(p::Union{AbstractArray, Tuple}, name)
    vars = @variables $(name)[axes(p)...]
    first(vars)
end

function symbolize(::Number, name)
    parameter(name)
end

function symbolize(p::Union{Symbolics.Num, Symbolics.Arr{Symbolics.Num}}, name)
    p
end

function symbolize(params::NamedTuple)
    NamedTuple{keys(params)}(Tuple(symbolize(v, Symbol("params_$(k)")) for (k,v) in pairs(params)))
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
