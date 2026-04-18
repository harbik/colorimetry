# Claude Code instructions for colorimetry-plot

This crate is a member of the `colorimetry` workspace. The primary instructions are in the
root [`CLAUDE.md`](../CLAUDE.md). All conventions, CI checks, and release rules defined
there apply here as well.

## Version sync

`colorimetry-plot` must always be at the same version as `colorimetry` and `colorimetry-cli`.
When bumping the version, update:

- `[package] version` in this file (`plot/Cargo.toml`)
- `colorimetry = {version = "..."}` dependency in this file (`plot/Cargo.toml`)
- The corresponding fields in `../Cargo.toml` and `../cli/Cargo.toml`

See the **Multi-crate version sync** section in the root `CLAUDE.md` for the full checklist,
including publication order and xtask commands to run before publishing.
