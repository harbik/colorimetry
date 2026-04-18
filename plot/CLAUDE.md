# Claude Code instructions for colorimetry-plot

This crate is a member of the `colorimetry` workspace. The primary instructions are in the
root [`CLAUDE.md`](../CLAUDE.md). All conventions, CI checks, and release rules defined
there apply here as well.

## Version sync

`colorimetry-plot` must always be at the same version as `colorimetry` and `colorimetry-cli`.
Version is managed via Cargo workspace inheritance — this crate uses `version.workspace = true`
and `colorimetry = {workspace = true, ...}`, so there is nothing to update here when bumping.

To bump the version, edit only `[workspace.package] version` and `[workspace.dependencies]
colorimetry version` in the root `../Cargo.toml`. See the **Multi-crate version sync** section
there for the full checklist, including publication order and xtask commands to run before
publishing.
