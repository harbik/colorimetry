// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2025, Harbers Bik LLC

//! # Rendable Trait Module
//!
//! This module defines the [`Rendable`] trait, which provides a unified interface for types that can be rendered as SVG elements with configurable view parameters.
//!
//! Types implementing [`Rendable`] can be composed into SVG documents, have their viewports adjusted, and be rendered to SVG output. The trait supports both required and optional methods, allowing for flexible and extensible rendering logic.
//!
//! ## Features
//! - Abstracts rendering logic for SVG-compatible types
//! - Supports view parameter configuration and adjustment
//! - Enables composition of complex SVG documents from multiple renderable components
//! - Provides optional methods for fine-grained control over view box position and size

use svg::node::element::SVG;

use crate::view::ViewParameters;

/// A trait for types that can be rendered as SVG elements with configurable view parameters.
///
/// Types implementing `Rendable` can be composed into SVG documents, have their viewports adjusted,
/// and be rendered to SVG output. This trait provides both required and optional methods for
/// interacting with the object's view parameters and rendering logic.
///
/// # Required Methods
/// - [`view_parameters`](Self::view_parameters): Returns the current view parameters for the object.
/// - [`set_view_parameters`](Self::set_view_parameters): Sets the view parameters for the object.
/// - [`render`](Self::render): Renders the object as an SVG element.
///
/// # Optional Methods
/// - [`set_x`](Self::set_x): Sets the x position of the view box.
/// - [`set_y`](Self::set_y): Sets the y position of the view box.
/// - [`width`](Self::width): Returns the width of the view box.
/// - [`height`](Self::height): Returns the height of the view box.
///
/// Implementors can override the optional methods for more efficient or specialized behavior.
pub trait Rendable {
    /// Returns the current view parameters for the object.
    fn view_parameters(&self) -> ViewParameters;

    /// Sets the view parameters for the object.
    fn set_view_parameters(&mut self, view_box: ViewParameters);

    /// Renders the object as an SVG element.
    fn render(&self) -> SVG;

    /// Sets the x position of the view box (optional).
    fn set_x(&mut self, x: i32) {
        let mut vb = self.view_parameters();
        vb.set_x(x);
        self.set_view_parameters(vb);
    }

    /// Sets the y position of the view box (optional).
    fn set_y(&mut self, y: i32) {
        let mut vb = self.view_parameters();
        vb.set_y(y);
        self.set_view_parameters(vb);
    }

    /// Returns the width of the view box (optional).
    fn width(&self) -> u32 {
        self.view_parameters().width()
    }

    /// Sets the width of the view box (optional).
    fn set_width(&mut self, width: u32) {
        let mut vb = self.view_parameters();
        vb.set_width(width);
        self.set_view_parameters(vb);
    }

    /// Returns the height of the view box (optional).
    fn height(&self) -> u32 {
        self.view_parameters().height()
    }

    /// Sets the height of the view box (optional).
    fn set_height(&mut self, height: u32) {
        let mut vb = self.view_parameters();
        vb.set_height(height);
        self.set_view_parameters(vb);
    }

    fn view_box(&mut self) -> (i32, i32, u32, u32) {
        let vp = self.view_parameters();
        vp.view_box()
    }

    fn set_view_box(&mut self, vx: i32, vy: i32, vw: u32, vh: u32) {
        let mut vp = self.view_parameters();
        vp.set_view_box(vx, vy, vw, vh);
        self.set_view_parameters(vp);
    }
}
