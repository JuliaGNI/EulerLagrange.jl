


function hamiltonian_variables(dimension::Int)
    t = parameter(:t)
    q = @variables (q(t))[1:dimension]
    p = @variables (p(t))[1:dimension]
    return (t,q,p)
end
