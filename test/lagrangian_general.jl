using EulerLagrange
using Symbolics
using Test

lvars = lagrangian_variables(2)

@variables t
@variables (x(t))[1:2]
@variables (v(t))[1:2]

@test isequal(t, lvars[1])
@test all([isequal(x[i], lvars[2][i]) for i in eachindex(x, lvars[2])])
@test all([isequal(v[i], lvars[3][i]) for i in eachindex(v, lvars[3])])
