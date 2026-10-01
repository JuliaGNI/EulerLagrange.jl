# Known issues

What is known to be broken or incomplete and is not fixed yet. Delete an entry when its issue is
fixed; the fix goes in `CHANGELOG.md`.

## Upstream

### K1 · The Documentation job fails because the docs environment does not resolve

- location: `docs/Project.toml:12`
- evidence: from the repository root, on Julia 1.11.9,
  `julia --startup-file=no --project=docs -e 'import Pkg; Pkg.resolve()'` fails with
  `Unsatisfiable requirements detected for package GeometricIntegrators [dcce2d33]`. In the General
  registry, every GeometricIntegrators release from 0.18 on bounds `GeometricBase = "0.14.8 - 0.14"`,
  and the package requires `GeometricBase = "0.15.0"`. The job heals when GeometricIntegrators
  registers a release for GeometricBase 0.15. Then a `workflow_dispatch` of `Documenter.yml` on
  `main` must be green, and a later pull request deletes this entry.
- kind: upstream
- found: 2026-10-02
