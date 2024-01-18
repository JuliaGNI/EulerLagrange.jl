using EulerLagrange
using EulerLagrange: NullParameters
using Symbolics
using Test


params  = (a = 1, b = [2,2], c = (3,3))
sparams = symbolize(params)

@test isequal(sparams, symbolize(sparams))

@variables a, b[1:2], c[1:2]

@test isequal(a, symbolize(a, :a))
@test isequal(b, symbolize(b, :b))
@test isequal(c, symbolize(c, :c))

@test isequal(sparams.a, a)
@test isequal(sparams.b, b)
@test isequal(sparams.c, c)

@test symbolize(nothing) == NamedTuple()
@test symbolize(NullParameters()) == NamedTuple()
