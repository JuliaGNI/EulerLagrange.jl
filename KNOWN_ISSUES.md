# Known issues

What is known to be broken or incomplete and is not fixed yet. Delete an entry when its issue is
fixed; the fix goes in `CHANGELOG.md`.

## Package

### K2 · Every field of `LagrangianSystem`, `DegenerateLagrangianSystem` and `HamiltonianSystem` is `Any`

- location: `src/lagrangian.jl:31`
- evidence: the structs at `src/lagrangian.jl:31`, `src/lagrangian_degenerate.jl:27` and
  `src/hamiltonian.jl:73` (at `723656d`) declare their 21 fields `L`/`H`, `t`, `x`/`q`, `v`/`p`,
  `parameters`, `equations` and `functions` with no type, or with `::Any`. The inner constructors
  store the symbolic Lagrangian or Hamiltonian (a `Num`), the symbolic time and coordinates, the
  parameter `NamedTuple`, a `NamedTuple` of symbolic expressions, and a `NamedTuple` of
  `RuntimeGeneratedFunction`s from `generate_code`. Every field access, and so every accessor
  (`lagrangian`, `hamiltonian`, `variables`, `parameters`, `equations`, `functions`), infers as
  `Any`, and a call through `functions(sys)` dispatches at run time. On Julia 1.13.3,
  `julia --startup-file=no concrete-fields.jl <package>` prints one line per field, such as
  `EulerLagrange.LagrangianSystem.L::Any`, and ends with `Verdict: FAIL — 21 violations`. The fix
  gives each struct one type parameter per stored field, so that `new` fixes the concrete types;
  it changes the types' signatures and is a refactor of the three structs.
- kind: defect
- found: 2026-11-09
