use std::{cell::RefCell, collections::HashMap};

use svg::{
    node::element::{ClipPath, Path, Style, Symbol},
    Document, Node,
};

const DEFAULT_CSS: &str = "
    .default {fill: lightgray; stroke: lightgray; stroke-width:1;}
    text.default  {color: lightgray; stroke: none; stroke-width: 0; font-size: 12pt; font-family: sans-serif;}
";

pub const NORTH : i32 = 90;
pub const SOUTH : i32 = 270;
pub const WEST: i32 = 180;
pub const EAST: i32 = 0;
pub const NORTH_WEST: i32 = 135;
pub const NORTH_EAST: i32 = 45;
pub const SOUTH_WEST: i32 = 225;
pub const SOUTH_EAST: i32 = 315;

use crate::layer::Layer;

/// Represents an SVG document with a specified width and height.  It contains clip paths, styles,
/// layers, and symbols that can be used to create SVG graphics.  Although you can directly
/// manipulate the SVG document, it is recommended to use higher level objects in this library such
/// as `Chart` struct for creating scaled charts with x and y axis.
pub struct SvgDocument {
    width: u32,
    height: u32, // [width, height]
    pub(super) clip_paths: RefCell<Vec<ClipPath>>,
    pub(super) styles: RefCell<HashMap<String, String>>, // selector, and elements
    pub(super) layers: RefCell<Vec<Layer>>,
    pub(super) symbols: RefCell<Vec<Symbol>>,
}

impl SvgDocument {
    /// Creates a new `SvgDocument` with the specified width and height.
    /// # Arguments
    /// * `width` - The width of the SVG document.
    /// * `height` - The height of the SVG document.
    /// # Returns
    /// A new instance of `SvgDocument`.
    pub fn new(width: u32, height: u32) -> Self {
        SvgDocument {
            width,
            height,
            clip_paths: RefCell::new(Vec::new()),
            styles: RefCell::new(HashMap::new()),
            layers: RefCell::new(Vec::new()),
            symbols: RefCell::new(Vec::new()),
        }
    }

    /// Adds a path to the SVG document as a clip path with the specified ID.
    pub fn add_clip_path(&self, id: String, path: &Path) {
        let clip = ClipPath::new().set("id", id).add(path.clone());
        let mut clips = self.clip_paths.borrow_mut();
        clips.push(clip);
    }

    /// Adds a symbol to the SVG document.
    pub fn add_symbol(&self, symbol: Symbol) {
        let mut symbols = self.symbols.borrow_mut();
        symbols.push(symbol);
    }

    pub fn add_css_rule(self, select: &str, style: &str) -> Self {
        {
            let mut styles = self.styles.borrow_mut();
            styles.insert(select.to_string(), style.to_string());
        }
        self
    }

    pub fn add_layer(&self, layer: Layer) {
        let mut layers = self.layers.borrow_mut();
        layers.push(layer);
    }

    /// Creates an SVG document with the specified width and height.
    pub fn svg(&self) -> Document {
        let mut doc = Document::new()
            .set("viewBox", (0, 0, self.width, self.height))
            .set("width", "100%")
            .set("height", "100%");

        let styles = self.styles.borrow_mut();
        let mut content = styles
            .iter()
            .map(|(selector, style)| format!("{} {{{}}}", selector, style))
            .collect::<Vec<String>>()
            .join("\n");
        content.push_str(DEFAULT_CSS);
        if !content.is_empty() {
            doc = doc.add(Style::new(content));
        }

        //  add definitions for clip paths and symbols
        let mut defs = svg::node::element::Definitions::new();

        let clips = self.clip_paths.borrow();
        for clip_path in clips.iter() {
            defs.append(clip_path.clone());
        }

        let symbols = self.symbols.borrow();
        for symbol in symbols.iter() {
            defs.append(symbol.clone());
        }

        doc.append(defs);

        let layers = self.layers.borrow();
        for layer in layers.iter() {
            doc = doc.add(layer.0.clone());
        }

        doc
    }

    pub fn render(&self) -> String {
        self.svg().to_string()
    }

    pub fn save(&self, filename: &str) -> Result<(), std::io::Error> {
        let doc = self.svg();
        svg::save(filename, &doc)?;
        Ok(())
    }
}
