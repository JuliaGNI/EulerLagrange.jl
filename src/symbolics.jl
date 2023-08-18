
parameter(name::Symbol) = Num(Sym{Real}(name))


function substitute_parameters(code, params)
    if length(params) > 0
        paramstr = string(Tuple("params_" * string(k) * ", " for k in keys(params))...)
        paramstr = paramstr[begin:prevind(paramstr, first((findlast(",", paramstr))))]
        
        code_str = string(code)
        code_str = replace(code_str, paramstr => "params")
        code_str = replace(code_str, "params_" => "params.")

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
