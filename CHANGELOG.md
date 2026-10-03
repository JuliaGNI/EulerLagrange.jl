# Release Notes

All notable changes to EulerLagrange.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31. 28 versions were released before it, the most recent
`v0.5.1`, and none of them are written up here: the record of that history is `git log` and the
tags. It is named as a gap rather than reconstructed, because a changelog assembled after the
fact loses exactly the reasoning that makes it worth keeping.

## [Unreleased] — targeting 0.5.3

### Changed

- The `[compat]` floor of RuntimeGeneratedFunctions is raised from `0.5` to `0.5.9`. Compat only:
  every Symbolics 7 release requires RuntimeGeneratedFunctions 0.5.9 or later, so 0.5.0 to 0.5.8
  could never install beside `Symbolics = "7"`, and the declared floor was false. The resolved
  environments do not change.

- The `symbolize(p, name)` method that returns a symbolic `p` unchanged has no type parameter
  any more. Its `where T` was unbound for a `Symbolics.Arr` argument, which Aqua reported (issue
  #27). The method accepts the same arguments and dispatch is unchanged. Aqua's unbound-argument
  check is no longer marked broken: Julia 1.14 does not report this method, so the broken mark
  failed the nightly job as an unexpected pass.

## [0.5.2]

### Changed

- The `[compat]` floors are raised to `GeometricBase = "0.15.0"`, `GeometricEquations = "0.21.5"`
  and `julia = "1.11"`, because GeometricBase 0.15 declares its stubs public and requires Julia 1.11.

- `test/lagrangian_solar_system.jl` no longer allocates the unread second momentum buffer `ṗ₂`,
  which clears fatou's one `unused-binding` finding in `test/`. The test checks the same values.

- Every tracked file is now Unicode NFC-normalised. Nine stored `ẋ` (43 times), `ṗ` (21), `ż` (15),
  `ḡ` (12), `ū` (7) and `ṽ` (6) as a base letter plus a combining mark, inherited from macOS rather
  than chosen. `q̇`, `v̄`, `f̄`, `f̃`, `p̃` and `ψ̃` have no precomposed codepoint and are unchanged.

  Nothing about the compiled code changes: Julia's parser normalises identifiers to NFC, so the
  symbols were already precomposed and dispatch, field names and method resolution are untouched.
  No string literal was affected, and the repository contains no doctests at all. What changes is
  that the source now matches what a keyboard, an editor search, a `grep` pattern or an automated
  replacement produces — in an NFD file a pattern typed in NFC matches nothing at all, silently.

  Every changed file is exactly the NFC normalisation of its predecessor, `docs/src/caveats.md`
  among them. Two further `ż` occurrences carry a tilde as well: `ż̃` composes its dot and keeps
  the tilde, there being no fully precomposed form, so a grep for `ż` in the old files found 17
  where 15 are counted above.

- Test suite reorganised. Dependencies moved from `[extras]`/`[targets]` to `test/Project.toml`
  (with `[sources]` for the package). Two wrapper files removed; seven test files now listed directly
  in `test/runtests.jl` with `@safetestset`. Total: ten existing test files (unchanged) plus new
  `test/quality/aqua.jl`. Aqua code-quality checks marked broken: issue #26 (ambiguities), #27
  (unbound type parameter). `Project.toml` gains `LinearAlgebra = "1"` compat entry, enforced by
  Aqua's deps_compat check (issue #28). `test/symbolize_tests.jl` renamed to `test/symbolics.jl`.
  No source changes.

- `test/Project.toml` no longer carries a `[compat]` entry for `GeometricEquations`, which the
  root `Project.toml` also depends on. The test environment contains the package, so the root's
  bound (`GeometricEquations = "0.21"`) already applies there; the removed `"0.21.4"` could only
  narrow it, and the tests then ran on narrower bounds than the package claims. The rule: the test
  and docs environments carry no `[compat]` entry for a dependency of the root. Test-only bounds
  are unchanged.
