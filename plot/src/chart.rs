use std::ops::{Deref, DerefMut};



use svg::{node::element::{path::Data, Group, Path}, Node};

use crate::{canvas::Canvas, transforms::TransformMatrix};


#[derive(Clone)]
pub struct Chart<'a> {
    pub canvas: &'a Canvas,
    pub group: Group,
    pub id: String,
    pub scale: [[f64; 2]; 2],
}

impl<'a> Chart<'a> {
    pub fn new(
        canvas: &'a Canvas,
        id: &str,
        target: [u32; 4],
        scale: [[f64; 2]; 2],
        clip: bool,
        style: Option<&str>,
    ) -> Self {
        let mut group = Group::new()
            .set("id", id)
            .set("vector-effect", "non-scaling-stroke")
            .set(
                "transform",
                TransformMatrix::new(target, scale).to_chart_string(),
            );
        let path = to_path(
            &format!("clip-{}", id),
            [
                (scale[0][0], scale[1][0]),
                (scale[0][1], scale[1][0]),
                (scale[0][1], scale[1][1]),
                (scale[0][0], scale[1][1]),
            ],
            true,
        );
        if clip {
            canvas.add_clip_path(format!("clip-{}", id), &path);
            group = group.set("clip-path", format!("url(#clip-{})", id));
        }
        if let Some(style) = style {
            let mut styles = canvas.styles.borrow_mut();
            styles.insert(format!("#{}", id), style.to_string());
            group = group.add(path)
                .set("id", format!("{}", id));
        }
        Chart {
            canvas,
            group,
            id: id.to_string(),
            scale,
        }
    }

    pub fn render(&self) {
        self.canvas.add_chart(self.clone());
    }

    pub fn draw_grid(&mut self, x_step: f64, y_step: f64, style: &str) {
        let mut data = Data::new();
        let [[xmin, xmax], [ymin, ymax]] = self.scale;
        let mut x = xmin;
        while x < xmax {
            data = data.move_to((x, ymin)).line_to((x, ymax));
            x += x_step;
        }
        let mut y = ymin;
        while y < ymax {
            data = data.move_to((xmin, y)).line_to((xmax, y));
            y += y_step;
        }
        let path = Path::new()
            .set("d", data)
            .set("vector-effect", "non-scaling-stroke")
            .set("style", style);
        self.group.append(path);
    }

    pub fn rectangle(&mut self, id: &str, x: f64, y: f64, width: f64, height: f64, style: String) {
        let rect = svg::node::element::Rectangle::new()
            .set("id", id)
            .set("x", x)
            .set("y", y)
            .set("width", width)
            .set("height", height);
        self.group.append(rect);
        let mut styles = self.canvas.styles.borrow_mut();
        styles.insert(format!("#{}", id), style);
    }

    pub fn path(
        &mut self,
        id: &str,
        data: impl IntoIterator<Item = (f64, f64)>,
        close: bool,
        style: Option<&str>,
    ) {
        let path = to_path(id, data, close);
        self.group.append(path);
        if let Some(style) = style {
            let mut styles = self.canvas.styles.borrow_mut();
            styles.insert(format!("#{}", id), style.to_string());
        }
    }

    pub fn add_style(&self, select: &str, style: &str) {
        let mut styles = self.canvas.styles.borrow_mut();
        styles.insert(select.to_string(), style.to_string());
    }
}

impl<'a> Deref for Chart<'a> {
    type Target = Group;

    fn deref(&self) -> &Self::Target {
        &self.group
    }
}

impl<'a> DerefMut for Chart<'a> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.group
    }
}

fn to_path(id: &str, data: impl IntoIterator<Item = (f64, f64)>, close: bool) -> Path {
    let mut path_data = Data::new();
    for xy in data {
        if path_data.is_empty() {
            path_data = path_data.move_to(xy);
        } else {
            path_data = path_data.line_to(xy);
        }
    }
    if close {
        path_data = path_data.close();
    }
    Path::new()
        .set("id", id.to_string())
        .set("d", path_data.clone())
        .set("vector-effect", "non-scaling-stroke")
}