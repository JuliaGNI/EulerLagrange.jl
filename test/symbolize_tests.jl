using EulerLagrange
using EulerLagrange: NullParameters
using Symbolics
using Test


params  = (a = 1, b = [2,2], c = (3,3))
sparams = symbolize(params)

@test isequal(sparams, symbolize(sparams))

@variables params_a, params_b[1:2], params_c[1:2]

@test isequal(params_a, symbolize(params_a, :params_a))
@test isequal(params_b, symbolize(params_b, :params_b))
@test isequal(params_c, symbolize(params_c, :params_c))

@test isequal(sparams.a, params_a)
@test isequal(sparams.b, params_b)
@test isequal(sparams.c, params_c)

@test symbolize(nothing) == NamedTuple()
@test symbolize(NullParameters()) == NamedTuple()
