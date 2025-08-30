
# Colorimetry Plots

<!-- cargo-rdme start -->

## Colorimetry Plot Library

This colormimetry library provides functionality for generating SVG-based color plots,
with both low-level and higher level API, based on top of the [Rust-SVG](https://crates.io/crates/svg) library.
It includes generating basic 2D (x,y) charts composed of several layers, and more complex chromaticity diagrams with
the spectral locus, and gamut fills.
Here is an example of a XY Chromaticity Diagram for the DisplayP3 color space, using the CIE 2015 observer:
![XY Chromaticity Diagram](https://harbik.github.io/colorimetry/img/srgb_gamut.svg)
Plots are built up in `Layers` using coordinate transformations between plot space and world coordinates.

### Modules

- `chart`: Chart composition and rendering.
- `layer`: Layered rendering support for compositing multiple plot elements.
- `rendable`: Traits and types for objects that can be rendered to SVG.
- `spectrum`: Spectrum data visualization and color representation.
- `style_attr`: SVG styling attributes and utilities.
- `svgdoc`: SVG document creation and manipulation.
- `view`: Viewport and plot area management.

### Core Features

- **Layered Architecture**: Build complex plots by compositing multiple layers
- **Coordinate Transforms**: Seamless conversion between plot and world coordinate systems
- **SVG Generation**: High-quality vector graphics output
- **Configurable Styling**: Flexible styling system for plot elements
- **Precision Control**: Configurable floating-point precision for clean output

<!-- cargo-rdme end -->
