
parameter(name::Symbol) = Num(Sym{Real}(name))
# variable(name::Symbol) = Num(Variable{Symbolics.FnType{Tuple{Any}, Real}}(name))
# variable(name::Symbol) = Num(variable(name; T = Symbolics.FnType{Tuple{Any}, Real}))
