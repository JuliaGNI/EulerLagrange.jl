# Caveats

The Hamiltonian or Lagrangian should be provided in a form as simple as possible to aid the symbolic algebra system.
Many expressions lead to involved symbolic expressions that can be hard to differentiate and lead to errors in the code generation pipeline.

For example, scalar products should be specified as `x ⋅ x` and not as `x' * x`.
As the transposition also entails complex conjugation, the resulting symbolic expression is substantially more complicated, which can cause problems down the line.

Another potential problem is an array operation, like the dot product, that is preceded by scalar multiplication like in `α * x ⋅ x`.
This can cause problems due to the order of operations, which is `(α * x) ⋅ x`, where first `x` is rescaled by `α`, implying the creation of a temporary array holding the result, which is then contracted with `x`.
This should rather be expressed as `x ⋅ x * α` or `α * (x ⋅ x)`, where the scalar product is computed first and then the result, which now is just a scalar, is multiplied by `α`.
