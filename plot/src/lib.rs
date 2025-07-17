
//! # Colorimetry Plot Library
//!
//! This module provides functionality for generating and manipulating SVG-based plots,
//! including axes, charts, chromaticity diagrams, layers, renderable objects, spectra,
//! SVG document handling, coordinate transforms, and view management.
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
pub mod axis;
pub mod chart;
pub mod chromaticity;
pub mod layer;
pub mod rendable;
pub mod spectrum;
pub mod svgdoc;
pub mod transforms;
pub mod view;

use std::sync::atomic::{AtomicUsize, Ordering};

static COUNTER: AtomicUsize = AtomicUsize::new(0);
static PRECISION: i32 = 1; // 1 decimal place

pub fn new_id() -> String {
    format!("id{}", COUNTER.fetch_add(1, Ordering::Relaxed))
}

pub fn last_id() -> String {
    COUNTER.load(Ordering::Relaxed).to_string()
}

pub fn assign_class_and_style<T: svg::Node>(
    node: &mut T,
    class: Option<&str>,
    style: Option<&str>,
) {
    if let Some(class) = class {
        node.assign("class", class);
    } else {
        node.assign("class", "default");
    }
    if let Some(style) = style {
        node.assign("style", style);
    }
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
