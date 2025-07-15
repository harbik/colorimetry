use std::{cell::RefCell, collections::HashMap};

use svg::{
    node::element::{ClipPath, Path, Style, Symbol, SVG},
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


/// Represents an SVG document with a specified width and height.  It contains clip paths, styles,
/// layers, and symbols that can be used to create SVG graphics.  Although you can directly
/// manipulate the SVG document, it is recommended to use higher level objects in this library such
/// as `Chart` struct for creating scaled charts with x and y axis.
pub struct SvgDocument {
    width: u32,
    height: u32, // [width, height]
    pub(super) clip_paths: RefCell<Vec<ClipPath>>,
    pub(super) styles: RefCell<HashMap<String, String>>, // selector, and elements
    pub(super) nodes: RefCell<Vec<Box<dyn Node>>>, // layers and use
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
            nodes: RefCell::new(Vec::new()), // layers, use, and svg elements
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
    pub fn add_symbol(&self, symbol: impl Into<Symbol>) {
        let mut symbols = self.symbols.borrow_mut();
        symbols.push(symbol.into());
    }

    pub fn add_css_rule(self, select: &str, style: &str) -> Self {
        {
            let mut styles = self.styles.borrow_mut();
            styles.insert(select.to_string(), style.to_string());
        }
        self
    }

    // TODO: check node type to select where to add it to
    /// Adds a node to the SVG document. The node can be any type that implements the `Node` trait.
    pub fn add(&self, node: impl Into<Box<dyn Node>>) {
        let a_node = node.into();
        let mut nodes = self.nodes.borrow_mut();
        nodes.push(a_node);
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

        let nodes = self.nodes.borrow();
        for node in nodes.iter() {
            doc = doc.add(node.clone());
        }

        doc
    }


    pub fn place(&mut self, svg: impl Into<SVG>, x: u32, y: u32, width: u32, height: u32) {
        self.add(svg.into()
            .set("x", x)
            .set("y", y)
            .set("width", width)
            .set("height", height)
            .set("preserveAspectRatio", "xMinYMin meet"));
    }
    
    /*
    pub fn render(&self) -> String {
        self.svg().to_string()
    }
     */

    pub fn save(&self, filename: &str) -> Result<(), std::io::Error> {
        let doc = self.svg();
        svg::save(filename, &doc)?;
        Ok(())
    }
}
