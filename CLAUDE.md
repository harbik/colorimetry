# Claude Code instructions for colorimetry

## Pre-release / CI checks

Whenever code is modified or tests are run, always execute the full xtask pipeline in this order and report any errors before marking work complete:

```
cargo xtask check   # fmt, clippy, build, rdme
cargo xtask test    # test --all-features and --no-default-features
cargo xtask doc     # rustdoc with --deny warnings
```

These three commands together are the definition of "the build is green" for this project.

## Coding conventions

- Prefer editing existing files over creating new ones.
- Do not add docstrings, comments, or type annotations to code that was not changed.
- Private helper functions (e.g. `rf_from_de`, `N_ANGLE_BIN`) must not be linked
  from public doc comments — rustdoc rejects them with `--deny warnings`. Use plain
  backtick notation instead: `` `rf_from_de` ``.
- Feature-gated code must compile both with and without the feature. The xtask test
  command runs both `--all-features` and `--no-default-features`.

## Reference standards

CIE and IES standard documents are stored locally in `standards/` (gitignored — never committed to
the repo). The PDFs are purchased standards that cannot be redistributed; each developer must
obtain their own copy. When implementing or verifying anything against a standard, read the
relevant document from that directory. Key documents to look for:

- `standards/Colorimetry, 4th Edition (CIE15-2018).pdf` — CIE 15:2018 Colorimetry
  (available for purchase at cie.co.at)

Add new documents to `standards/` and list them here (with purchase URL) as they are acquired.

## Domain context

- The library implements CIE and ANSI/IES colorimetry standards.
- The colorimetric calculations should be implemented in accordance with the CIE 15 standard, in particular CIE 15:2004 and CIE 15:2018.
- CFI (Color Fidelity Index) follows **CIE 224:2017** / **ANSI/IES TM-30-20/24**.
  - 99 CES samples, 16 hue bins, CIECAM02-UCS J'a'b' colour space.
  - Bin assignment always uses the **reference-source** hue angle, not the test-source.
  - Per-bin metrics: Rf,hj (fidelity), Rcs,hj (chroma shift fraction), Rhs,hj (hue shift, rad).
- Tolerances for per-bin tests are intentionally wide (±10 Rf, ±0.05 Rcs, ±0.04 rad)
  because a single sample near a bin boundary can shift a bin centroid significantly.
- Reference data for tests: IES TM-30 Spectral Calculator (ies.org) for F12;
  LuxPy (`spd_to_ies_tm30_metrics`) for F1/F2 as a baseline only.

## Release process

Before opening or merging a PR, the following must all pass with no errors or warnings:

1. `cargo xtask check`   — fmt, clippy (-D warnings), build, rdme
2. `cargo xtask test`    — all-features and no-default-features
3. `cargo xtask doc`     — rustdoc with --deny warnings

The CHANGELOG.md should be updated for any user-facing change.
The version in Cargo.toml follows semver; the library is pre-1.0 so breaking changes
are allowed but must be noted in the changelog.
Deprecated items use `#[deprecated(since = "x.y.z", note = "use X instead")]`.

## Multi-crate version sync

This workspace contains three published crates that must always share the same version number:

| Crate | Cargo.toml |
|---|---|
| `colorimetry` | `Cargo.toml` |
| `colorimetry-plot` | `plot/Cargo.toml` |
| `colorimetry-cli` | `cli/Cargo.toml` |

When bumping the version for a release, update **all** of the following in a single commit:

1. `Cargo.toml` — `[package] version`
2. `plot/Cargo.toml` — `[package] version` **and** the `colorimetry = {version = "..."}` dependency
3. `cli/Cargo.toml` — `[package] version` **and** the `colorimetry = {version = "..."}` dependency
4. Any install instructions in `README.md`, `cli/README.md`, or `plot/README.md` that pin a specific version

Before publishing, run the full xtask pipeline from the workspace root (which covers all member crates):

```sh
cargo xtask check
cargo xtask test
cargo xtask doc
```

Publish in dependency order: `colorimetry` first, then `colorimetry-plot` and `colorimetry-cli`
(both depend on `colorimetry`, so the new version must be on crates.io before the dependents can be published).

## Tests

- Use `assert!(approx::abs_diff_eq!(got, want, epsilon = X), "bin {j}: got {got}, want {want}")`
  rather than `assert_abs_diff_eq!` with a message argument (the macro does not accept one).
- Tests that depend on optional features must be gated with `#[cfg(feature = "...")]`.
- Do not mock internals. Tests hit real computation paths.
