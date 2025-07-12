use std::{cell::RefCell, collections::HashMap};

use svg::{
    node::element::{ClipPath, Path, Style, Symbol},
    Document, Node,
};

use crate::layer::Layer;


pub struct SvgDocument {
    width: u32,
    height: u32, // [width, height]
    pub(super) clip_paths: RefCell<Vec<ClipPath>>,
    pub(super) styles: RefCell<HashMap<String, String>>, // selector, and elements
    pub(super) layers: RefCell<Vec<Layer>>,
    pub(super) symbols: RefCell<Vec<Symbol>>
}

impl SvgDocument {
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
    
    /*
    pub fn render_chart(&self, chart: Chart) {
        self.add_group(chart.xyplot.clone());
        self.add_group(chart.annotations.clone());
        self.add_group(chart.axes.clone());
    }
     */

    pub fn add_clip_path(&self, id: String, path: &Path) {
        let clip = ClipPath::new().set("id", id).add(path.clone());
        let mut clips = self.clip_paths.borrow_mut();
        clips.push(clip);
    }

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
        let mut layers= self.layers.borrow_mut();
        layers.push(layer);
    }

    /// Creates an SVG document with the specified width and height.
    pub fn svg(&self) -> Document {
        let mut doc = Document::new()
            .set("viewBox", (0, 0, self.width, self.height))
            .set("width", "100%")
            .set("height", "100%");

        // add styles
        let styles = self.styles.borrow();
        let content = styles
            .iter()
            .map(|(selector, style)| format!("{} {{{}}}", selector, style))
            .collect::<Vec<String>>()
            .join("\n");
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

