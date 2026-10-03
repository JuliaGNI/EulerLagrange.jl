using Aqua
using EulerLagrange
using Test

Aqua.test_all(EulerLagrange;
    ambiguities = (broken = true,))     # issue #26
