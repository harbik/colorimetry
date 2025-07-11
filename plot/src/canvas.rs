use std::{cell::RefCell, collections::HashMap};

use svg::{
    node::element::{ClipPath, Group, Path, Style, Symbol},
    Document, Node,
};

use crate::{chart::{XYChart, to_path}, new_id, transforms::TransformMatrix};

pub struct Canvas {
    width: u32,
    height: u32, // [width, height]
    pub(super) clip_paths: RefCell<Vec<ClipPath>>,
    pub(super) styles: RefCell<HashMap<String, String>>, // selector, and elements
    pub(super) groups: RefCell<Vec<Group>>,
    pub(super) symbols: RefCell<Vec<Symbol>>
}

impl Canvas {
    pub fn new(width: u32, height: u32) -> Self {
        Canvas {
            width,
            height,
            clip_paths: RefCell::new(Vec::new()),
            styles: RefCell::new(HashMap::new()),
            groups: RefCell::new(Vec::new()),
            symbols: RefCell::new(Vec::new()),
        }
    }
    
    /// Creates a new chart on the canvas.
    ///
    /// The chart is positioned using a target rectangle defined as `[left, top, width, height]`,
    /// and scaled using a coordinate system defined as `[[x_min, x_max], [y_min, y_max]]`.
    ///
    /// An optional CSS `class` and `style` can also be applied.
    pub fn add_chart<'a>(
        &'a self,
        target: [u32; 4], // left, top, width, height
        scale: [[f64; 2]; 2],
        class: Option<&'a str>,
        style: Option<&'a str>,
    ) -> XYChart<'a> {
        // chart alwyas gets an id
        let id = new_id();
        let transform_matrix = TransformMatrix::new(target, scale);
        let mut group = Group::new().set(
            "transform",
            transform_matrix.to_chart_string(),
        );
        let mut path = to_path(
            [
                (scale[0][0], scale[1][0]),
                (scale[0][1], scale[1][0]),
                (scale[0][1], scale[1][1]),
                (scale[0][0], scale[1][1]),
            ],
            true,
        );
        // chart always gets a clip path
        self.add_clip_path(format!("clip-{}", id), &path);
        group = group.set("clip-path", format!("url(#clip-{})", id));
        // don't add a background if no style or class is given
        if style.is_some() || class.is_some() {
            if let Some(style) = style {
                 path = path.set("style", style);
            }
            if let Some(class) = class {
                path = path.set("class", class);
            }
            group.append(path);
        }
        let annotations = Group::new()
            .set("id", format!("annotations-{}", id))
            .set("class", "annotations");

        let axes = Group::new()
            .set("id", format!("axes-{}", id))
            .set("class", "axes");
        XYChart {
            canvas: &self,
            xyplot: group,
            annotations,
            axes,
            scale,
            target,
            transform_matrix

        }
    }


    pub fn render_chart(&self, chart: XYChart) {
        self.add_group(chart.xyplot.clone());
        self.add_group(chart.annotations.clone());
        self.add_group(chart.axes.clone());
    }

    pub fn add_clip_path(&self, id: String, path: &Path) {
        let clip = ClipPath::new().set("id", id).add(path.clone());
        let mut clips = self.clip_paths.borrow_mut();
        clips.push(clip);
    }

    pub fn add_symbol(&self, symbol: Symbol) {
        let mut symbols = self.symbols.borrow_mut();
        symbols.push(symbol);
    }

    pub fn add_style(self, select: &str, style: &str) -> Self {
        {
            let mut styles = self.styles.borrow_mut();
            styles.insert(select.to_string(), style.to_string());
        }
        self
    }

    pub fn add_group(&self, group: Group) {
        let mut groups = self.groups.borrow_mut();
        groups.push(group);
    }

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

        let groups = self.groups.borrow();
        for group in groups.iter() {
            doc = doc.add(group.clone());
        }

        doc
    }

    pub fn save(&self, filename: &str) -> Result<(), std::io::Error> {
        let doc = self.svg();
        svg::save(filename, &doc)?;
        Ok(())
    }
}

