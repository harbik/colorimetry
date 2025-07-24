//! # Colorimetry Plot Library
//!
//! This colormimetry library provides functionality for generating SVG-based color plots,
//! with both low-level and higher lelvel API, based on top of the Rust-SVG library.
//! It includes generating basic 2D (x,y) charts composed of several layers, and more complex chromaticity diagrams with
//! the spectral locus, and gamut fills.
//! Plots are build up in `Layers` using transformations `to_plot` and `to_world`.
//!
//! ## Modules
//!
//! - `axis`: Axis rendering and management.
//! - `chart`: Chart composition and rendering.
//! - `chromaticity`: Chromaticity diagram utilities.
//! - `layer`: Layered rendering support.
//! - `rendable`: Traits and types for renderable objects.
//! - `spectrum`: Spectrum data and visualization.
//! - `svgdoc`: SVG document creation and manipulation.
//! - `transforms`: Coordinate and geometric transforms.
//! - `view`: Viewport and plot size management.
//!
//! ## Utilities
//!
//! - Unique ID generation for SVG elements.
//! - Class and style assignment for SVG nodes.
//! - Floating-point rounding utilities with configurable precision.
//!
//! ## Usage
//!
//! Import the desired modules and use the provided functions to construct and manipulate SVG plots.
//!
//! ## License
//!
//! This library is dual-licensed under the MIT License and the Apache License (Version 2.0).
//! You may choose either license when using this library.
pub mod chart;
pub mod layer;
pub mod rendable;
pub mod spectrum;
pub mod svgdoc;
pub mod view;
pub mod style_attr;

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
