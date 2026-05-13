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
- `standards/CIE 2017 Colour Fidelity Index for accurate scientific use (CIE224-2017).pdf` — CIE 224:2017
  Colour Fidelity Index for accurate scientific use (available for purchase at cie.co.at)
- `standards/TM-30-20-E1.pdf` - ANSI/IES TM-30-20 Technical Memorandum: IES Method for evaluating light source color rendition (available for purchase at <https://store.ies.org/product-category/science>).

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

This is the complete checklist for cutting a new release. All published artifacts must carry the
same version number: the three Rust crates on crates.io, the `colorimetry` npm package (WASM), and
the GitHub release tag.

### 1. Update CHANGELOG.md

Move every entry under `## [Unreleased]` into a new dated section, e.g.:

```markdown
## [0.0.9] - 2026-04-20
```

Add a diff link at the bottom of the file following the existing pattern. All three changelogs
(`CHANGELOG.md`, `cli/CHANGELOG.md`, `plot/CHANGELOG.md`) should be updated if they contain
unreleased entries.

### 2. Bump version numbers

All of the following must change from the old version to the new version in **one commit**:

| File | What to change |
|---|---|
| `Cargo.toml` | `[workspace.package] version` |
| `Cargo.toml` | `[workspace.dependencies] colorimetry` version pin |
| `pkg/package.json` | `"version"` field (npm / WASM package) |
| `README.md` | any `colorimetry = "x.y.z"` install snippets |
| `cli/README.md` | any version-pinned install snippets |
| `plot/README.md` | any version-pinned install snippets |

The three Rust crates (`colorimetry`, `colorimetry-plot`, `colorimetry-cli`) inherit version via
`version.workspace = true` — only `Cargo.toml` at the workspace root needs to change.

### 3. Run the full xtask pipeline

```sh
cargo xtask check   # fmt, clippy (-D warnings), build, rdme
cargo xtask test    # --all-features and --no-default-features
cargo xtask doc     # rustdoc with --deny warnings
```

All three must pass with zero errors and zero warnings before proceeding.

### 4. Rebuild the WASM package

```sh
cargo xtask wasm    # regenerates pkg/ via wasm-pack + wasm-opt
```

Commit any changes to `pkg/` (`.wasm`, `.js`, `.d.ts` files) together with the version bump commit,
or as a follow-up commit before tagging.

### 5. Commit and tag

```sh
git add -p          # stage version bumps, CHANGELOG, pkg/ changes
git commit -m "chore: release v0.0.9"
git tag v0.0.9
git push origin main --tags
```

The tag triggers the GitHub release. Create a GitHub Release from the tag (via the web UI or
`gh release create v0.0.9 --notes-from-tag`) and paste the relevant CHANGELOG section as the
release notes.

### 6. Publish to crates.io

`colorimetry-plot` depends on the external `cmx` crate (at `../cmx`), which in turn depends on
`colorimetry`. Because `cmx` pins a specific `colorimetry` version, it must be updated and
published **before** `colorimetry-plot` can be published.

#### 6a. Update and publish `cmx` (at `../cmx`)

```sh
# In ../cmx:
# 1. Bump colorimetry version in Cargo.toml to the new version
# 2. Fix any compilation errors caused by renamed APIs
# 3. Run tests: cargo test
# 4. Commit and publish:
git add Cargo.toml Cargo.lock <any changed src files>
git commit -m "chore: bump colorimetry to x.y.z, release vA.B.C"
cargo publish
```

#### 6b. Update `colorimetry-plot` to the new `cmx` version

```sh
# In plot/Cargo.toml, bump cmx to the version just published
# Then commit:
git add plot/Cargo.toml Cargo.lock
git commit -m "chore: update colorimetry-plot to use cmx vA.B.C"
git push origin main
```

#### 6c. Publish all three crates in dependency order

```sh
cargo publish -p colorimetry
# wait for crates.io to index it (usually ~30 seconds), then:
cargo publish -p colorimetry-plot
cargo publish -p colorimetry-cli
```

### 7. Publish to npm

```sh
cd pkg
npm publish --access public
cd ..
```

Verify the new version appears on npmjs.com before closing the release.

### Per-PR changelog and deprecation notes

- `CHANGELOG.md` should be updated for any user-facing change, not just at release time.
- The version in `Cargo.toml` follows semver; the library is pre-1.0 so breaking changes are
  allowed but must be noted in the changelog.
- Deprecated items use `#[deprecated(since = "x.y.z", note = "use X instead")]`.

## Tests

- Use `assert!(approx::abs_diff_eq!(got, want, epsilon = X), "bin {j}: got {got}, want {want}")`
  rather than `assert_abs_diff_eq!` with a message argument (the macro does not accept one).
- Tests that depend on optional features must be gated with `#[cfg(feature = "...")]`.
- Do not mock internals. Tests hit real computation paths.
