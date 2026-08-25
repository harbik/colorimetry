# Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## Categories each change fall into

* **Added**: for new features.
* **Changed**: for changes in existing functionality.
* **Deprecated**: for soon-to-be removed features.
* **Removed**: for now removed features.
* **Fixed**: for any bug fixes.
* **Security**: in case of vulnerabilities.

## Unreleased

## [0.1.0] - 2026-08-25

### Changed

* Released in lockstep with `colorimetry` 0.1.0 and built against it. The library release
  carries breaking changes (the `spectral-io` feature flag was removed and `MeasurementKind`
  became `MeasurementType`) — see the [library changelog](../CHANGELOG.md). No command,
  argument, or output of the `color` binary changed.

## [0.0.9] - 2026-04-20

### Changed

* The version is now inherited from the workspace (`version.workspace = true`), so the CLI,
  the library, and the WASM package always release under one number. This release exists to
  keep that alignment; it carries no user-facing CLI changes.
* Internal: adapted to the library's `XYZ::values()` -> `XYZ::to_array()` rename. Command
  output and arguments are unaffected.

## [0.0.8] - 2025-08-28

### Added

* Introduced `clap`-based structures and methods
* Added some basic commands, mostly to test the framework
* In this version, only `color sample` is implemented
