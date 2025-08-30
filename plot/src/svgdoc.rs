// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2025, Harbers Bik LLC

//! SVG document representation and rendering.
//!
//! This module provides the `SvgDocument` struct, which allows you to create and manipulate SVG documents
//! and render them to SVG files. It includes methods for adding clip paths, symbols, and plots,
//! as well as setting the document's width, height, and margin.

use std::process;

use svg::{
    node::element::{ClipPath, Group, Path, Rectangle, Style, Symbol},
    Document, Node,
};

use crate::{rendable::Rendable, view::ViewParameters};

/// Default CSS styles for the SVG document.
/// This CSS is applied to all elements in the SVG document by default.
/// It sets the default fill, stroke, stroke-width, and font styles for text elements.
/// You can append additional styles using the `append_scss` method of the `SvgDocument`.
/// The styles are compiled to CSS using the `grass` crate when the SVG document is rendered.
const DEFAULT_CSS: &str = "
    * { 
        fill:none; 
        stroke:black; 
        stroke-width:1px;
        stroke-linecap: round;
    }
    text  {
        fill:black; 
        stroke: none; 
        font-size: 16px; 
        font-family: Helvetica, Arial, sans-serif; 
    }
";

pub const NORTH: i32 = 90;
pub const SOUTH: i32 = 270;
pub const WEST: i32 = 180;
pub const EAST: i32 = 0;
pub const NORTH_WEST: i32 = 135;
pub const NORTH_EAST: i32 = 45;
pub const SOUTH_WEST: i32 = 225;
pub const SOUTH_EAST: i32 = 315;

/// Represents an SVG document with a specified width and height.  It contains clip paths, styles,
/// layers, and symbols that can be used to create SVG graphics.  Although you can directly
/// manipulate the SVG document, it is recommended to use higher level objects in this library such
/// as `Chart` struct for creating scaled charts with x and y axis.
///
/// The `SvgDocument` struct is the main entry point for creating SVG documents in this library.
///
/// A SVG-global style sheet can be set using the `scss` field, which is a string containing SCSS
/// styles. The styles are compiled to CSS using the `grass` crate when the SVG document is rendered.
/// This allows for flexible styling of SVG elements using SCSS syntax, and verified the validity of the SCSS parameters
/// before rendering. You can append additional SCSS styles using the `append_scss` method.
///
/// Any css file is also a valid SCSS file, so you can use any CSS file as a stylesheet.
///
/// The style sheet is embedded in the SVG document as a `<style>` element, which is
/// the most convenient way to handle SVG images with styles.
///
/// You don't have to use a global style, you can set styles for individual elements using the
/// `StyleAttr` parameter in all the plot element and plot functions.
#[derive(Default)]
pub struct SvgDocument {
    pub(super) scss: String, // SCSS stylesheet content
    pub(super) clip_paths: Vec<ClipPath>,
    pub(super) symbols: Vec<Symbol>,
    pub(super) plots: Vec<Box<dyn Rendable>>,
    pub(super) nodes: Vec<Box<dyn Node>>, // layers and use
    pub(super) margin: u32,
    pub(super) width: Option<u32>,
    pub(super) height: Option<u32>,
}

impl SvgDocument {
    const DEFAULT_MARGIN: u32 = 50;
    /// Creates a new `SvgDocument` with default settings only,
    /// without any nodes, symbols, or plots.
    pub fn new() -> Self {
        SvgDocument {
            clip_paths: Vec::new(),
            scss: DEFAULT_CSS.to_string(),
            nodes: Vec::new(), // layers, use, and svg elements
            symbols: Vec::new(),
            plots: Vec::new(),
            margin: Self::DEFAULT_MARGIN,
            width: None,
            height: None,
        }
    }

    /// Appends SCSS styles to the document's stylesheet.
    /// This method allows you to add additional styles to the existing SCSS content.
    /// The styles are compiled to CSS when the document is rendered.
    ///
    /// Most covenient is to use a separate SCSS file and include it using the `append_scss` method
    /// using `include_str!` macro, which reads the file at compile time.
    /// The SCSS styles are compiled to CSS using the `grass` crate when the SVG document is rendered,
    /// and embedded in the resulting SVG document as a `<style>` element.
    ///
    /// Using a seperate SCSS file allows you to keep your styles organized and maintainable,
    /// and also allows you to edit the styles in a lanaguate aware editor with syntax highlighting and autocompletion.
    ///
    /// # Example
    /// ```ignore
    /// // The file in the include_str! macro is just an example, and gives a compile error if the file does not exist.
    /// use colorimetry_plot::svgdoc::SvgDocument;
    ///
    /// // Create a new SVG document and append SCSS styles from a file
    /// let svg_doc = SvgDocument::new()
    ///     .append_scss(include_str!("path/to/your/styles.scss"));
    /// ```
    pub fn append_scss(mut self, css: &str) -> Self {
        self.scss.push_str(css);
        self
    }

    /// Sets the margin for the SVG document.
    /// The margin is applied to the width and height of the document when calculating the size of sub plots.
    /// The default margin is 50 pixels.
    pub fn set_margin(mut self, margin: u32) -> Self {
        self.margin = margin;
        self
    }

    /// Sets the width of the SVG document.
    /// If not set, the width will be calculated based on the size of the sub plots and the margin.
    /// If the width is larger than the maximum allowed size, an error will be printed and the program will exit.
    /// The width is set in pixels.
    pub fn set_width(mut self, width: u32) -> Self {
        self.width = Some(width);
        self
    }

    /// Sets the height of the SVG document.
    /// If not set, the height will be calculated based on the size of the sub plots and the margin.
    /// If the height is larger than the maximum allowed size, an error will be printed and the program will exit.
    /// The height is set in pixels.
    pub fn set_height(mut self, height: u32) -> Self {
        self.height = Some(height);
        self
    }

    /// Adds a path to the SVG document as a clip path with the specified ID.
    pub fn add_clip_path(mut self, id: String, path: &Path) -> Self {
        let clip = ClipPath::new().set("id", id).add(path.clone());
        self.clip_paths.push(clip);
        self
    }

    /// Adds a symbol to the SVG document.
    pub fn add_symbol(mut self, symbol: impl Into<Symbol>) -> Self {
        self.symbols.push(symbol.into());
        self
    }

    /// Adds a node to the SVG document, typically used for headers, footers, or other SVG elements.
    /// They will be added to the SVG document as-is, without any transformations or scaling,
    /// and placed in front of any other content, into a layer with the id "front".
    pub fn add_node(mut self, node: impl Into<Box<dyn Node>>) -> Self {
        self.nodes.push(node.into());
        self
    }

    /// Adds a sub plot to the SVG document.
    /// The sub plot must implement the `Rendable` trait, which allows it to be rendered as an SVG element.
    /// The sub plot will be positioned based on the document's flow, which is calculated by the `flow` method.
    /// If the sub plot is larger than the document size, an error will be printed and the program will exit.
    pub fn add_plot(mut self, svg_sub: Box<dyn Rendable>) -> Self {
        self.plots.push(svg_sub);
        self
    }

    /// Calculates the positions of the sub plots on the document.
    /// The positions are calculated based on the document's width and height,
    /// and the size of the sub plots. If there is only one sub plot, it will be centered in the document.
    /// If there are multiple sub plots, the flow is not yet implemented and will return `todo!()`.
    /// The sub plots will not be scaled, but positioned based on the document's flow.
    /// If the sub plot is larger than the document size, an error will be printed and the program will exit.
    /// # Returns
    /// A vector of tuples containing the x and y positions of each sub plot in the document.
    pub fn flow(&self) -> Vec<(u32, u32)> {
        match self.plots.len() {
            1 => {
                let svg_sub = &self.plots[0];
                let doc_width = self.width();
                let doc_height = self.height();

                let width = svg_sub.width();
                let height = svg_sub.height();

                if width > doc_width || height > doc_height {
                    eprintln!("Error: SVG sub-element is larger than the document size.");
                    eprintln!(
                        "Document size: {doc_width}x{doc_height}, Sub-element size: {width}x{height}");
                    eprintln!("Please adjust the size of the SVG sub-element or the document.");
                    process::exit(1); // Exit with error code 1
                }

                let x = doc_width / 2 - width / 2;
                let y = doc_height / 2 - height / 2;
                vec![(x, y)]
            }
            _ => todo!(),
        }
    }

    pub fn calculate_subplots_size_with_margin(&self) -> (u32, u32) {
        match self.plots.len() {
            1 => {
                let svg_sub = &self.plots[0];
                (
                    svg_sub.width() + 2 * self.margin,
                    svg_sub.height() + 2 * self.margin,
                )
            }
            _ => (800, 600), // Default size for the SVG document
        }
    }

    pub fn save(self, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
        Ok(svg::save(filename, &self.render())?)
    }
}

impl Rendable for SvgDocument {
    fn view_parameters(&self) -> ViewParameters {
        let (subs_width, subs_height) = self.calculate_subplots_size_with_margin();
        let width = self.width.unwrap_or(subs_width);
        let height = self.height.unwrap_or(subs_height);
        ViewParameters::new(0, 0, width, height, width, height)
    }

    fn set_view_parameters(&mut self, _view_box: ViewParameters) {
        /* do nothing, calculated up on rendering */
    }

    fn render(&self) -> Document {
        let vp = self.view_parameters();
        let mut doc = Document::new()
            .set("viewBox", vp.to_string())
            .set("width", vp.width())
            .set("height", vp.height())
            .set("xmlns", "http://www.w3.org/2000/svg")
            .set("xmlns:xlink", "http://www.w3.org/1999/xlink")
            .set("version", "1.1")
            .set("class", "colorimetry-plot");

        //  let scss_content = format!("{}\n{}", DEFAULT_CSS, self.scss);
        let css_content = match grass::from_string(self.scss.clone(), &grass::Options::default()) {
            Ok(css) => css,
            Err(e) => {
                eprintln!("Failed to compile SCSS: {e}");
                process::exit(1); // Exit with error code 1
            }
        };
        doc = doc.add(Style::new(css_content));

        //  add definitions for clip paths and symbols
        let mut defs = svg::node::element::Definitions::new();

        for clip_path in self.clip_paths.iter() {
            defs.append(clip_path.clone());
        }

        self.symbols.iter().for_each(|symbol| {
            defs.append(symbol.clone());
        });
        doc.append(defs);

        let background = Group::new()
            .set("id", "background")
            .set("class", "background")
            .add(
                Rectangle::new()
                    .set("x", 0)
                    .set("y", 0)
                    .set("width", vp.width())
                    .set("height", vp.height()),
            );
        doc = doc.add(background);

        // add plots
        for (plot, (x, y)) in self.plots.iter().zip(self.flow()) {
            let rendered_plot = plot.render().set("x", x).set("y", y);
            doc = doc.add(rendered_plot);
        }

        // add all items on the svg document
        for node in self.nodes.iter() {
            doc = doc.add(node.clone());
        }

        doc
    }
}
