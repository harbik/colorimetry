
# Colorimetry Command Line Interface Tool

<!-- cargo-rdme start -->

This `colorimetry-cli` crate provides the core functionality for the `color` command-line tool, enabling users to convert color data, create plots, analyze color spaces, and perform other color-related tasks directly from the terminal.

## Installation

Ensure you have Rust and Cargo installed. You can install them using [rustup](https://rustup.rs).

To install `colorimetry-cli`, run:

```bash
cargo install colorimetry-cli
```

To verify your installation and see available commands and options, run:

```bash
color --help
```

## Usage Example

To convert the spectral data in the `sample1.csv` file to CIELAB values, use the `color` subcommand:

```bash
color sample sample1.csv -m lab
```

The tool currently supports CSV files with wavelengths in the first column and spectral values in the second column.  
Use the `-m` option to specify the target color model, such as `xyz`, `srgb`, or `lab`.

> **Note:** In code, use `colorimetry_cli` as the crate name, since Rust replaces dashes with underscores.

<!-- cargo-rdme end -->
