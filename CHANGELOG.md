# Release Notes

All notable changes to EulerLagrange.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31 and deliberately holds no entries. 28 versions were
released before it, the most recent `v0.5.1`, and none of them are written up here: the
record of that history is `git log` and the tags. It is named as a gap rather than
reconstructed, because a changelog assembled after the fact loses exactly the reasoning that
makes it worth keeping. The `[Unreleased]` target below is provisional — confirm it when the
first entry is written.

## [Unreleased] — targeting 0.6.0

### Changed

- Every tracked source file is now Unicode NFC-normalised. Nine files stored `ẋ` (43 times), `ṗ`
  (21), `ż` (15), `ḡ` (12), `ū` (7) and `ṽ` (6) as a base letter plus a combining mark, inherited
  from macOS rather than chosen.

  Nothing about the compiled code changes: Julia's parser normalises identifiers to NFC, so the
  symbols were already precomposed and dispatch, field names and method resolution are untouched.
  No string literal was affected, and no changed line falls inside a doctest block. What changes is
  that the source now matches what a keyboard, an editor search, a `grep` pattern or an automated
  replacement produces — in an NFD file a pattern typed in NFC matches nothing at all, silently.

  Every changed file is exactly the NFC normalisation of its predecessor, `docs/src/caveats.md`
  among them. `ż̃` keeps its combining tilde, having no fully precomposed form.

### New Features

### Bug Fixes

### Breaking Changes

## Open Issues
