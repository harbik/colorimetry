use std::{collections::HashMap};

use svg::{
    node::element::{ClipPath, Path, Style, Symbol},
    Document, Node,
};

use crate::{view:: ViewParameters, rendable::Rendable};

const DEFAULT_CSS: &str = "
    .default {fill: lightgray; stroke: lightgray; stroke-width:1;}
    text.default  {fill: lightgray; stroke: none; stroke-width: 0; font-size: 16px; font-family: san-serif;}
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
#[derive(Default)]
pub struct SvgDocument {
    pub(super) view_parameters: ViewParameters,
    pub(super) styles: HashMap<String, String>, // selector, and elements
    pub(super) clip_paths: Vec<ClipPath>,
    pub(super) symbols: Vec<Symbol>,
    pub(super) plots: Vec<Box<dyn Rendable>>,
    pub(super) nodes: Vec<Box<dyn Node>>, // layers and use
}

impl SvgDocument {
    /// Creates a new `SvgDocument` with the specified width and height.
    /// # Arguments
    /// * `width` - The width of the SVG document.
    /// * `height` - The height of the SVG document.
    /// # Returns
    /// A new instance of `SvgDocument`.
    pub fn new(width: u32, height: u32) -> Self {
        let mut view_parameters = ViewParameters::default();
        view_parameters.set_width(width);
        view_parameters.set_height(height);
        view_parameters.set_view_box(0, 0, width, height);

        SvgDocument {
            view_parameters,
            clip_paths: Vec::new(),
            styles: HashMap::new(),
            nodes: Vec::new(), // layers, use, and svg elements
            symbols: Vec::new(),
            plots: Vec::new(),

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

/*

    /// Creates an SVG document with the specified width and height.

    pub fn place(&mut self, svg_sub: impl Rendable) {
        self.view_parameters.extend(&svg_sub.view_parameters());
        self.add(svg_sub.render());
    }

    pub fn place_at_position(&mut self, svg_sub: impl Rendable, x: i32, y: i32) {
        self.view_parameters.extend_with_pos(&svg_sub.view_parameters(), x, y);
        let sub_svg = svg_sub.render()
            .set("x", x)
            .set("y", y);
        self.add(sub_svg);
    }
 */

    pub fn add_svg(&mut self, svg_sub: Box<dyn Rendable>) {
        self.plots.push(svg_sub);
    }

    /*

    // TODO TODO Need to set height and width of the subplot!!!
    // 
    pub fn place_center(&mut self, svg_sub: impl Rendable) {
        let width = svg_sub.width();
        let height = svg_sub.height();
        let doc_width = self.width();
        let doc_height = self.height();
        let x = doc_width /2 - width /2;
        let y = doc_height /2 - height /2;

        let sub_svg = svg_sub.render()
            .set("x", x)
            .set("y", y);
        self.add(sub_svg);
    }
     */

    pub fn position(&self) -> Vec<(u32, u32)> {
        match self.plots.len() {
            1 =>  {
                let svg_sub = &self.plots[0];
                let doc_width = self.width();
                let doc_height = self.height();

                let width = svg_sub.width();
                let height = svg_sub.height();
                let x = doc_width /2 - width /2;
                let y = doc_height /2 - height /2;
                vec!((x, y))
            }
            _ => todo!()
        }
    }

    pub fn save(&mut self, filename: &str) -> Result<(), std::io::Error> {
        /*
        match self.plots.len() {
            0 => {
                // If there are no plots, just render the document
                return svg::save(filename, &self.render());
            }
            1 => {
                let doc_width = self.width();
                let doc_height = self.height();

                let svg_sub = &mut self.plots[0];
                let width = svg_sub.width();
                let height = svg_sub.height();
                let x = doc_width /2 - width /2;
                let y = doc_height /2 - height /2;
                svg_sub.set_x(x as i32);
                svg_sub.set_y(y as i32);

                /*
                let sub_svg_rendered = svg_sub.render()
                    .set("x", x)
                    .set("y", y);

                let mut doc = self.render();
                doc = doc.add(sub_svg_rendered);
                 */

         */
        return svg::save(filename, &self.render());
    }
}

impl Rendable for SvgDocument {
    fn view_parameters(&self) -> ViewParameters {
        self.view_parameters.clone()
    }

    fn set_view_parameters(&mut self, view_box: ViewParameters) {
        self.view_parameters = view_box;
    }
    
    fn render(&self) -> Document {
        let mut doc = Document::new()
            .set("viewBox", self.view_parameters.to_string())
            .set("width", self.view_parameters.width())
            .set("height", self.view_parameters.height());

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

        // add plots
        for (plot, (x, y)) in self.plots.iter().zip(self.position()) {
            let rendered_plot = 
            plot
                .render()
                .set("x", x)
                .set("y", y);
            doc = doc.add(rendered_plot);
        }

        for node in self.nodes.iter() {
            doc = doc.add(node.clone());
        }

        doc
    }
}