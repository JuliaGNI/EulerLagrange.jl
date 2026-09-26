using Aqua
using EulerLagrange
using Test

Aqua.test_all(EulerLagrange;
    ambiguities = (broken = true,),     # issue #26
    unbound_args = (broken = true,),    # issue #27
    deps_compat = (broken = true,))     # issue #28, no compat entry for LinearAlgebra
