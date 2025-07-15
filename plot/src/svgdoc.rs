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
    margin: u32,
    pub(super) placements: Vec<(i32, i32, u32, u32)>, // (x, y, width, height)
    pub(super) clip_paths: Vec<ClipPath>,
    pub(super) styles: HashMap<String, String>, // selector, and elements
    pub(super) nodes: Vec<Box<dyn Node>>, // layers and use
    pub(super) symbols: Vec<Symbol>,
}

impl SvgDocument {
    /// Creates a new `SvgDocument` with the specified width and height.
    /// # Arguments
    /// * `width` - The width of the SVG document.
    /// * `height` - The height of the SVG document.
    /// # Returns
    /// A new instance of `SvgDocument`.
    pub fn new(margin: u32) -> Self {
        SvgDocument {
            margin,
            placements: Vec::new(),
            clip_paths: Vec::new(),
            styles: HashMap::new(),
            nodes: Vec::new(), // layers, use, and svg elements
            symbols: Vec::new(),
        }
    }

    /// Adds a path to the SVG document as a clip path with the specified ID.
    pub fn add_clip_path(&mut self, id: String, path: &Path) {
        let clip = ClipPath::new().set("id", id).add(path.clone());
        self.clip_paths.push(clip);
    }

    /// Adds a symbol to the SVG document.
    pub fn add_symbol(&mut self, symbol: impl Into<Symbol>) {
        self.symbols.push(symbol.into());
    }

    pub fn add_css_rule(mut self, select: &str, style: &str) -> Self{
        self.styles.insert(select.to_string(), style.to_string());
        self
    }

    // TODO: check node type to select where to add it to
    /// Adds a node to the SVG document. The node can be any type that implements the `Node` trait.
    pub fn add(&mut self, node: impl Into<Box<dyn Node>>) {
        self.nodes.push(node.into());
    }


    /// Creates an SVG document with the specified width and height.
    pub fn svg(&self) -> Document {
        let mut doc = Document::new()
          //  .set("viewBox", (0, 0, self.width, self.height))
            .set("width", "100%")
            .set("height", "100%");

        let mut content = self.styles
            .iter()
            .map(|(selector, style)| format!("{} {{{}}}", selector, style))
            .collect::<Vec<String>>()
            .join("\n");

        content.push_str(DEFAULT_CSS);
        doc = doc.add(Style::new(content));

        //  add definitions for clip paths and symbols
        let mut defs = svg::node::element::Definitions::new();

        for clip_path in self.clip_paths.iter() {
            defs.append(clip_path.clone());
        }

        self.symbols.iter().for_each(|symbol| {
            defs.append(symbol.clone());
        });

        doc.append(defs);

        for node in self.nodes.iter() {
            doc = doc.add(node.clone());
        }

        doc
    }


    pub fn place(&mut self, svg: impl Into<SVG>, x: u32, y: u32, width: Option<u32>, height: Option<u32>) {
        let mut svg = svg.into();
        svg.assign("x",x);
        svg.assign("y",y);
        if let Some(w) = width {
            svg.assign("width", w);
        }
        if let Some(h) = height {
            svg.assign("heigth", h);
        }
        self.add(svg)
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
