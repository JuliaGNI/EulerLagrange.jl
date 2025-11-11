using EulerLagrange
using EulerLagrange: NullParameters
using Symbolics
using Test


params = (a=1, b=[2, 2], c=(3, 3))
sparams = symbolize(params)

@test isequal(sparams, symbolize(sparams))

@variables aₚ bₚ[1:2] cₚ[1:2]

@test isequal(aₚ, symbolize(aₚ, :aₚ))
@test isequal(bₚ, symbolize(bₚ, :bₚ))
@test isequal(cₚ, symbolize(cₚ, :cₚ))

@test isequal(sparams.a, aₚ)
@test isequal(sparams.b, bₚ)
@test isequal(sparams.c, cₚ)

@test symbolize(nothing) == NamedTuple()
@test symbolize(NullParameters()) == NamedTuple()
