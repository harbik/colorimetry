use std::{cell::RefCell, collections::HashMap};

use svg::{node::element::{ClipPath, Group, Path, Style}, Document};

use crate::chart::Chart;


pub struct Canvas {
    width: u32,
    height: u32, // [width, height]
    pub(super) clip_paths: RefCell<Vec<ClipPath>>,
    pub(super) styles: RefCell<HashMap<String, String>>, // selector, and elements
    pub(super) groups: RefCell<HashMap<String, Group>>,
}

impl Canvas {
    pub fn new(width: u32, height: u32) -> Self {
        Canvas {
            width,
            height,
            clip_paths: RefCell::new(Vec::new()),
            styles: RefCell::new(HashMap::new()),
            groups: RefCell::new(HashMap::new()),
        }
    }

    pub fn add_chart(&self, chart: Chart) {
        self.add_group(&chart.id, chart.group.clone());
    }

    pub fn add_clip_path(&self, id: String, path: &Path) {
        let clip = ClipPath::new().set("id", id).add(path.clone());
        let mut clips = self.clip_paths.borrow_mut();
        clips.push(clip);
    }

    pub fn add_style(self, select: &str, style: &str) {
        let mut styles = self.styles.borrow_mut();
        styles.insert(select.to_string(), style.to_string());
    }

    pub fn add_group(&self, id: &str, group: Group) {
        let mut groups = self.groups.borrow_mut();
        groups.insert(id.to_string(), group);
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

        let clips = self.clip_paths.borrow();
        for clip_path in clips.iter() {
            doc = doc.add(clip_path.clone());
        }

        let groups = self.groups.borrow();
        for (id, group) in groups.iter() {
            doc = doc.add(group.clone().set("id", id.clone()));
        }

        doc
    }

    pub fn save(&self, filename: &str) -> Result<(), std::io::Error> {
        let doc = self.svg();
        svg::save(filename, &doc)?;
        Ok(())
    }
}
