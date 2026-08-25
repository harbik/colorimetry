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
## [0.1.0] - 2026-08-25
```

Add a diff link at the bottom of the file following the existing pattern.

**Every released version gets a dated section in BOTH changelogs** (`CHANGELOG.md` and
`cli/CHANGELOG.md`) — not only the ones that happen to have accumulated `Unreleased`
entries. The CLI inherits its version from the workspace, so it ships on every release
whether or not its own code changed. When nothing user-facing changed, say exactly that:

```markdown
## [0.0.9] - 2026-04-20

### Changed

* Version-alignment release: no changes to the `color` binary's commands, arguments,
  or output.
```

"No user-facing change" is information a reader wants; a missing section just looks like
an oversight. In 0.1.0 both CLI releases (0.0.9 and 0.1.0) were reconstructed from git
history after the fact because an empty `## Unreleased` was read as "nothing to do".

Before writing the sections, confirm nothing user-facing is unrecorded:

```sh
git log <previous-tag>..HEAD --oneline            # library
git log <previous-tag>..HEAD --oneline -- cli/    # CLI
```

**Date the section with the day you actually publish, not the day you prepare it.** The
0.1.0 section was written as `2026-05-15` and shipped on `2026-08-25`; re-check the date
immediately before tagging (step 5) and correct it if the release slipped.

### 2. Bump version numbers

All of the following must change from the old version to the new version in **one commit**:

| File | What to change |
|---|---|
| `Cargo.toml` | `[workspace.package] version` |
| `Cargo.toml` | `[workspace.dependencies] colorimetry` version pin |
| `deno.jsonc` | `"version"` field **and** the `jsr:@harbik/colorimetry@^x.y.z` import pin (JSR package) |
| `README.md` | any `colorimetry = "x.y"` install snippets |
| `cli/README.md` | any version-pinned install snippets |

Both Rust crates (`colorimetry`, `colorimetry-cli`) inherit version via
`version.workspace = true` — only `Cargo.toml` at the workspace root needs to change.

`pkg/package.json` is **not** edited by hand and is git-ignored (`pkg/.gitignore` is `*`).
wasm-pack regenerates it from `Cargo.toml` on every build, so its version follows from
step 4 — verify it reads the new version with `npm pack --dry-run` before publishing.

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
or as a follow-up commit before tagging. Note that `colorimetry.d.ts` embeds the rustdoc comments
of every `#[wasm_bindgen]` item, so doc-only edits to `src/**/wasm.rs` also require this rebuild.

### 5. Commit and tag

```sh
git add -p          # stage version bumps, CHANGELOG, pkg/ changes
git commit -m "chore: release v0.1.0"
git tag v0.1.0
git push origin main --tags
```

The tag triggers the GitHub release. Create a GitHub Release from the tag (via the web UI or
`gh release create v0.1.0 --notes-from-tag`) and paste the relevant CHANGELOG section as the
release notes.

### 6. Publish to crates.io

Publish both workspace crates in dependency order:

```sh
cargo publish -p colorimetry
# wait for crates.io to index it (usually ~30 seconds), then:
cargo publish -p colorimetry-cli
```

**Note:** `colorimetry-plot` is now a standalone crate at
<https://github.com/harbik/colorimetry-plot>. It has its own release process. It depends
on `colorimetry` with a loose `"0.1"` pin, so patch releases of colorimetry do not
require a `colorimetry-plot` release.

### 7. Publish to npm

```sh
cd pkg
npm publish --access public
cd ..
```

Verify the new version appears on npmjs.com before closing the release.

### 8. Publish to JSR (Deno)

```sh
deno publish --dry-run
deno publish
```

Publishes `@harbik/colorimetry` from `deno.jsonc`, which exports the same `pkg/` build.
JSR versions are immutable — check <https://jsr.io/@harbik/colorimetry> first.

### Per-PR changelog and deprecation notes

- `CHANGELOG.md` should be updated for any user-facing change, not just at release time.
  A change that touches `cli/` goes in `cli/CHANGELOG.md` under its own `## Unreleased`.
  Doing this per-PR is what keeps step 1 of the release a formality instead of an
  archaeology exercise.
- The version in `Cargo.toml` follows semver; the library is pre-1.0 so breaking changes are
  allowed but must be noted in the changelog.
- Deprecated items use `#[deprecated(since = "x.y.z", note = "use X instead")]`.

## Tests

- Use `assert!(approx::abs_diff_eq!(got, want, epsilon = X), "bin {j}: got {got}, want {want}")`
  rather than `assert_abs_diff_eq!` with a message argument (the macro does not accept one).
- Tests that depend on optional features must be gated with `#[cfg(feature = "...")]`.
- Do not mock internals. Tests hit real computation paths.
