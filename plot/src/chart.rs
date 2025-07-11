use std::{fmt::Display, ops::{Deref, DerefMut}};

use svg::{
    node::{element::{path::Data, Group, Path, Text}},
    Node,
};

use crate::{axis::{Axis, AxisSide}, canvas::Canvas, last_id, transforms::TransformMatrix};

#[derive(Clone)]
pub struct XYChart<'a> {
    pub canvas: &'a Canvas,
    pub xyplot: Group,
    pub axes: Group,
    pub annotations: Group,
    pub scale: [[f64; 2]; 2],
    pub target: [u32; 4], // left, top, width, height
    pub transform_matrix: TransformMatrix,
}

impl<'a> XYChart<'a> {

    pub fn add_axis(mut self, side: AxisSide, step: f64, show_labels: bool, class: Option<&str>) -> Self {
        let min_max = match side {
            AxisSide::Bottom | AxisSide::Top => self.scale[0],
            AxisSide::Left | AxisSide::Right => self.scale[1],
        };
        let x_axis = Axis::new(self.target, min_max, step, side, show_labels, class);
        self.axes.append(Group::from(x_axis));
        self
    }

    pub fn add_x_axis(self, x_step: f64, margin: u32, class: Option<&str>, style: Option<&str>) -> Self {
        let [left, top, width, height] = self.target;
        let [[xmin, xmax], [_ymin, _ymax]] = self.scale;

        let mut x_axis = Group::new();
        let mut x = xmin;
        while x<xmax {
            let x_pos = left as f64 + (x - xmin) / (xmax - xmin) * width as f64;
            let y_pos = top as f64 + height as f64 + margin as f64;
            x_axis.append(
                Text::new(format!("{:.2}", x))
                .set("x", x_pos)
                .set("y", y_pos)
                .set("text-anchor", "middle")
                .set("dominant-baseline", "hanging")
            );
            x += x_step;
        }
        if let Some(class) = class {
            x_axis = x_axis.set("class", class);
        }
        if let Some(style) = style {
            x_axis = x_axis.set("style", style);
        }
        self.canvas.add_group(x_axis);
        self
    }

    pub fn add_y_axis(self, y_step: f64, margin: u32, class: Option<&str>, style: Option<&str>) -> Self {
        let [left, top, _width, height] = self.target;
        let [[_xmin, _xmax], [ymin, ymax]] = self.scale;

        let mut y_axis = Group::new();
        let mut y = ymin;
        while y < ymax {
            let txt = format!("{:.2}",  y );
            let txt_width = txt.chars().count() as f64 * 0.6 * 12.0; // Approximate width of text 
            let y_pos = (top + height) as f64 - (y - ymin) / (ymax - ymin) * height as f64;
            let x_pos = left as f64 - margin as f64 - txt_width as f64;
            y_axis.append(
                Text::new(txt)
                .set("x", x_pos)
                .set("y", y_pos)
                .set("text-anchor", "right")
                .set("dominant-baseline", "middle")
            );
            y += y_step;
        }
        if let Some(class) = class {
            y_axis = y_axis.set("class", class);
        }
        if let Some(style) = style {
            y_axis = y_axis.set("style", style);
        }
        self.canvas.add_group(y_axis);
        self
    }

    pub fn render(self) {
        self.canvas.render_chart(self);
    }

    pub fn draw_grid(self, x_step: f64, y_step: f64, class: Option<&str>, style: Option<&str>) -> Self {
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
        self.data(data, class, style)
    }

    // Add a style to the canvas, for the last created element without a style parameter.
    pub fn apply_style(self, style: &str) -> Self {
        let id = last_id();
        self.add_style(&format!("#{}", id), style)
    }

    pub fn draw_rect(
        mut self,
        id: &str,
        x: f64,
        y: f64,
        width: f64,
        height: f64,
        style: String,
    ) -> Self {
        let rect = svg::node::element::Rectangle::new()
            .set("id", id)
            .set("x", x)
            .set("y", y)
            .set("width", width)
            .set("height", height);
        self.xyplot.append(rect);
        let mut styles = self.canvas.styles.borrow_mut();
        styles.insert(format!("#{}", id), style);
        self
    }

    pub fn draw_dot(
        mut self,
        cx: f64,
        cy: f64,
        r: f64,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let (cx, cy) = self.transform_matrix.canvas(cx, cy);
        let mut circle = svg::node::element::Circle::new()
            .set("cx", cx)
            .set("cy", cy)
            .set("r", r);
        if let Some(class) = class {
            circle = circle.set("class", class);
        }
        if let Some(style) = style {
            circle = circle.set("style", style);
        }
        self.annotations.append(circle);
        self
    }

    pub fn draw_path(
        mut self,
        mut path: Path,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        if let Some(style) = style {
            path = path.set("style", style);
        }
        if let Some(class) = class {
            path = path.set("class", class);
        }
        self.xyplot.append(path);
        self
    }

    // Add Data, composed of e.g. move_to and line_to operations, to the Chart
    pub fn data(
        self,
        data: Data,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let path = Path::new()
            .set("d", data)
            .set("vector-effect", "non-scaling-stroke");
        self.draw_path(path, class, style)
    }

    /// Add a path to the chart, using coordinates from an iterator.
    pub fn draw_line(
        self,
        data: impl IntoIterator<Item = (f64, f64)>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        self.draw_path(to_path(data, false), class, style)
    }

    /// Add a path that connects the coordinates in an interator and closes the path.
    pub fn draw_area(
        self,
        data: impl IntoIterator<Item = (f64, f64)>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        self.draw_path(to_path(data, true), class, style)
    }

    pub fn add_style(self, select: &str, style: &str) -> Self {
        let mut styles = self.canvas.styles.borrow_mut();
        styles.insert(select.to_string(), style.to_string());
        self
    }
}

impl<'a> From<XYChart<'a>> for Group {
    fn from(chart: XYChart) -> Self {
        let group = Group::new()
            .set("id", format!("xyplot-{}", last_id()))
            .set("class", "xyplot")
            .add(chart.xyplot)
            .add(chart.annotations)
            .add(chart.axes);
       // if !chart.canvas.styles.borrow().is_empty() {
       //     group = group.set("style", chart.canvas.styles.borrow().to_string());
       // }
        group
    }
}

impl Display for XYChart<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.xyplot)
    }
}

impl<'a> Deref for XYChart<'a> {
    type Target = Group;

    fn deref(&self) -> &Self::Target {
        &self.xyplot
    }
}

impl<'a> DerefMut for XYChart<'a> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.xyplot
    }
}

pub(super) fn to_path(data: impl IntoIterator<Item = (f64, f64)>, close: bool) -> Path {
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
        //  .set("id", id.to_string())
        .set("d", path_data.clone())
        .set("vector-effect", "non-scaling-stroke")
}
