
# Colorimetry Plots

<!-- cargo-rdme start -->

## Colorimetry Plot Library

This colormimetry library provides functionality for generating SVG-based color plots,
with both low-level and higher lelvel API, based on top of the Rust-SVG library.
It includes generating basic 2D (x,y) charts composed of several layers, and more complex chromaticity diagrams with
the spectral locus, and gamut fills.
Plots are build up in `Layers` using transformations `to_plot` and `to_world`.

### Modules

- `axis`: Axis rendering and management.
- `chart`: Chart composition and rendering.
- `chromaticity`: Chromaticity diagram utilities.
- `layer`: Layered rendering support.
- `rendable`: Traits and types for renderable objects.
- `spectrum`: Spectrum data and visualization.
- `svgdoc`: SVG document creation and manipulation.
- `transforms`: Coordinate and geometric transforms.
- `view`: Viewport and plot size management.

### Utilities

- Unique ID generation for SVG elements.
- Class and style assignment for SVG nodes.
- Floating-point rounding utilities with configurable precision.

### Usage

Import the desired modules and use the provided functions to construct and manipulate SVG plots.

### License

This library is dual-licensed under the MIT License and the Apache License (Version 2.0).
You may choose either license when using this library.

<!-- cargo-rdme end -->
