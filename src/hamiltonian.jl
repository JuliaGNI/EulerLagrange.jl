


function hamiltonian_variables(dimension::Int)
    t = parameter(:t)
    
    @variables (q(t))[1:dimension]
    @variables (p(t))[1:dimension]
    
    return (t,q,p)
end
