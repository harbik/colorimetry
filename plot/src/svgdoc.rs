use std::process;

use svg::{
    node::element::{ClipPath, Path, Style, Symbol},
    Document, Node,
};

use crate::{rendable::Rendable, view::ViewParameters};

const DEFAULT_CSS: &str = "
    * { fill:none; stroke:black; stroke-width:1px; }
    text  {fill:black; stroke: none; stroke-width: 0; font-size: 16px; font-family: san-serif; }
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

    pub fn append_scss(mut self, css: &str) -> Self {
        self.scss.push_str(css);
        self
    }

    pub fn set_margin(mut self, margin: u32) -> Self {
        self.margin = margin;
        self
    }

    pub fn set_width(mut self, width: u32) -> Self {
        self.width = Some(width);
        self
    }

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

    /// Adds a node to the SVG document. The node can be any type that implements the `Node` trait.
    pub fn add_node(mut self, node: impl Into<Box<dyn Node>>) -> Self {
        self.nodes.push(node.into());
        self
    }

    pub fn add_svg(mut self, svg_sub: Box<dyn Rendable>) -> Self {
        self.plots.push(svg_sub);
        self
    }

    /// Calculates the positions of the sub plots on the document.
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

    pub fn subs_size_with_margin(&self) -> (u32, u32) {
        match self.plots.len() {
            1 => {
                let svg_sub = &self.plots[0];
                (svg_sub.width() + 2 * self.margin, svg_sub.height() + 2 * self.margin)
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
        let (subs_width, subs_height) = self.subs_size_with_margin();
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
