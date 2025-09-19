// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2025, Harbers Bik LLC

//! # Colorimetry Plot Library
//!
//! This colormimetry library provides functionality for generating SVG-based color plots,
//! with both low-level and higher level API, based on top of the [Rust-SVG](https://crates.io/crates/svg) library.
//! It includes generating basic 2D (x,y) charts composed of several layers, and more complex chromaticity diagrams with
//! the spectral locus, and gamut fills.
//! Here is an example of a XY Chromaticity Diagram for the DisplayP3 color space, using the CIE 2015 observer:
//! ![XY Chromaticity Diagram](https://harbik.github.io/colorimetry/img/displayp3_gamut.svg)
//! Plots are built up in `Layers` using coordinate transformations between plot space and world coordinates.
//!
//! ## Modules
//!
//! - `chart`: Chart composition and rendering.
//! - `layer`: Layered rendering support for compositing multiple plot elements.
//! - `rendable`: Traits and types for objects that can be rendered to SVG.
//! - `spectrum`: Spectrum data visualization and color representation.
//! - `style_attr`: SVG styling attributes and utilities.
//! - `svgdoc`: SVG document creation and manipulation.
//! - `view`: Viewport and plot area management.
//!
//! ## Core Features
//!
//! - **Layered Architecture**: Build complex plots by compositing multiple layers
//! - **Coordinate Transforms**: Seamless conversion between plot and world coordinate systems
//! - **SVG Generation**: High-quality vector graphics output
//! - **Configurable Styling**: Flexible styling system for plot elements
//! - **Precision Control**: Configurable floating-point precision for clean output

pub mod chart;
pub mod layer;
pub mod rendable;
pub mod style_attr;
pub mod svgdoc;
pub mod view;

pub use crate::style_attr::{class, id, style, StyleAttr};

use std::sync::atomic::{AtomicUsize, Ordering};

static COUNTER: AtomicUsize = AtomicUsize::new(0);
static PRECISION: i32 = 1; // 1 decimal place

/// Generates a new unique ID for SVG elements.
/// The ID is prefixed with "id" and is incremented each time this function is called.
/// This is useful for ensuring that each SVG element has a unique identifier.
/// # Returns
/// A unique ID string in the format "idN", where N is an incrementing number.
pub fn new_id() -> String {
    format!("id{}", COUNTER.fetch_add(1, Ordering::Relaxed))
}

/// Returns the last generated ID as a string.
/// This can be useful for referencing the last created SVG element without generating a new ID.
pub fn last_id() -> String {
    COUNTER.load(Ordering::Relaxed).to_string()
}

/// Rounds a floating-point value to the specified precision.
/// For example, round_to_precision(1.23456, 100.0) returns 1.23.
///
/// # Arguments
/// * `value` - The value to round.
/// * `precision` - The precision
///
/// # Returns
/// The rounded value.
pub fn round_to_precision(value: f64, precision: i32) -> f64 {
    let multiplier = 10f64.powi(precision);
    (value * multiplier).round() / multiplier
}

/// Rounds a floating-point value to the default precision.
/// For example, round_to_default_precision(1.23456) returns 1.2.
/// /// # Arguments
/// * `value` - The value to round.
///     
/// # Returns
/// The rounded value.
pub fn round_to_default_precision(value: f64) -> f64 {
    round_to_precision(value, PRECISION)
}
