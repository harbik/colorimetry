//! # Colorimetry Plot Library
//!
//! This colormimetry library provides functionality for generating SVG-based color plots,
//! with both low-level and higher level API, based on top of the Rust-SVG library.
//! It includes generating basic 2D (x,y) charts composed of several layers, and more complex chromaticity diagrams with
//! the spectral locus, and gamut fills.
//!
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
pub mod spectrum;
pub mod style_attr;
pub mod svgdoc;
pub mod view;

pub use crate::style_attr::StyleAttr;

use std::sync::atomic::{AtomicUsize, Ordering};

static COUNTER: AtomicUsize = AtomicUsize::new(0);
static PRECISION: i32 = 1; // 1 decimal place

pub fn new_id() -> String {
    format!("id{}", COUNTER.fetch_add(1, Ordering::Relaxed))
}

pub fn last_id() -> String {
    COUNTER.load(Ordering::Relaxed).to_string()
}

/// Sets the class and style attributes on an SVG node.
///
/// # Arguments
/// * `node` - The SVG node to modify.
/// * `class` - Optional class name to set. Defaults to "default" if None.
/// * `style` - Optional style string to set. No style is set if None.
///
/// # Returns
/// The modified SVG node.
pub fn set_class_and_style<T>(mut node: T, class: Option<&str>, style: Option<&str>) -> T
where
    T: svg::Node,
{
    node.assign("class", class.unwrap_or("default"));

    if let Some(style_value) = style {
        node.assign("style", style_value);
    }
    node
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
