use std::{
    fmt::Display,
    ops::{Deref, DerefMut},
};

use svg::{
    node::element::{path::Data, Group, Line, Path, Text},
    Node,
};

use crate::{
    axis::{Axis, AxisSide},
    last_id,
    layer::Layer,
    new_id,
    svgdoc::SvgDocument,
    transforms::TransformMatrix,
};

#[derive(Clone)]
pub struct Chart<'a> {
    pub svgdoc: &'a SvgDocument,
    pub xyplot: Layer,
    pub axes: Layer,
    pub annotations: Layer,
    pub scale: [[f64; 2]; 2],
    pub target: [u32; 4], // left, top, width, height
    pub transform_matrix: TransformMatrix,
}

impl<'a> Chart<'a> {
    pub const ANNOTATE_SEP: u32 = 2;


    /// Creates a new chart on the canvas.
    ///
    /// The chart is positioned using a target rectangle defined as `[left, top, width, height]`,
    /// and scaled using a coordinate system defined as `[[x_min, x_max], [y_min, y_max]]`.
    ///
    /// An optional CSS `class` and `style` can also be applied.
    pub fn new(
        svgdoc: &'a SvgDocument,
        target: [u32; 4], // left, top, width, height
        scale: [[f64; 2]; 2],
        class: Option<&'a str>,
        style: Option<&'a str>,
    ) -> Chart<'a> {
        // chart alwyas gets an id
        let id = new_id();
        let transform_matrix = TransformMatrix::new(target, scale);
        let mut scaled_layer = Layer::new().set("transform", transform_matrix.to_chart_string());
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
        svgdoc.add_clip_path(format!("clip-{}", id), &path);
        scaled_layer = scaled_layer.set("clip-path", format!("url(#clip-{})", id));
        // don't add a background if no style or class is given
        if style.is_some() || class.is_some() {
            if let Some(style) = style {
                path = path.set("style", style);
            }
            if let Some(class) = class {
                path = path.set("class", class);
            }
            scaled_layer.append(path);
        }
        let annotations = Layer::new()
            .set("id", format!("annotations-{}", id))
            .set("class", "annotations");

        let axes = Layer::new()
            .set("id", format!("axes-{}", id))
            .set("class", "axes");
        Chart {
            svgdoc,
            xyplot: scaled_layer,
            annotations,
            axes,
            scale,
            target,
            transform_matrix,
        }
    }

    pub fn add_axis(
        mut self,
        description: Option<&str>,
        side: AxisSide,
        step: f64,
        tick_length: u32,
        show_labels: bool,
        class: Option<&str>,
    ) -> Self {
        let min_max = match side {
            AxisSide::Bottom | AxisSide::Top => self.scale[0],
            AxisSide::Left | AxisSide::Right => self.scale[1],
        };
        let x_axis = Axis::new(
            description,
            self.target,
            min_max,
            step,
            side,
            tick_length,
            show_labels,
            class,
        );
        self.axes.append(Group::from(x_axis));
        self
    }

    pub fn render(&self) {
        let chart_layer = Layer::new()
            .set("id", format!("chart-{}", new_id()))
            .set("class", "chart")
            .add(self.axes.clone())
            .add(self.xyplot.clone())
            .add(self.annotations.clone());
        self.svgdoc.add_layer(chart_layer);
    }

    pub fn draw_grid(
        self,
        x_step: f64,
        y_step: f64,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
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

    pub fn draw_rect(
        mut self,
        id: &str,
        x: f64,
        y: f64,
        width: f64,
        height: f64,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let mut rect = svg::node::element::Rectangle::new()
            .set("id", id)
            .set("x", x)
            .set("y", y)
            .set("width", width)
            .set("height", height);

        if let Some(class) = class {
            rect.assign("class", class);
        }
        if let Some(style) = style {
            rect.assign("style", style);
        }
        self.xyplot.append(rect);
        self
    }

    pub fn annotate(
        mut self,
        cxy : (f64, f64),
        r: f64,
        len: u32,
        angle: i32,
        text: impl AsRef<str>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let angle = 360 - ((angle + 360) % 360) as i32;
        let (cx, cy) = cxy;
        let (cx, cy) = self.transform_matrix.canvas(cx, cy);
        let dx = len as f64 * (angle as f64).to_radians().cos();
        let dy = len as f64 * (angle as f64).to_radians().sin();
        let circle = svg::node::element::Circle::new()
            .set("cx", cx)
            .set("cy", cy)
            .set("r", r);
        let line = Line::new()
            .set("x1", cx)
            .set("y1", cy)
            .set("x2", cx as f64 + dx)
            .set("y2", cy as f64 + dy);

        let dxt = (len + Self::ANNOTATE_SEP) as f64 * (angle as f64).to_radians().cos();
        let dyt = (len + Self::ANNOTATE_SEP) as f64 * (angle as f64).to_radians().sin();
        let mut text = Text::new(text.as_ref())
            .set("x", cx as f64 + dxt)
            .set("y", cy as f64 + dyt);

        if let Some(class) = class {
            text.assign("class", class);
        } else {
            text.assign("class", "default");
        }

        match angle {
            0..55 | 305..=360 => { // easat
                text.assign("text-anchor", "start");
                text.assign("dominant-baseline", "middle");
            }
            55..125 => { // north
                text.assign("text-anchor", "middle");
                text.assign("dominant-baseline", "text-before-edge");
            }
            125..235 => { // west
                text.assign("text-anchor", "end");
                text.assign("dominant-baseline", "middle");
            }
            _ => { // south
                text.assign("text-anchor", "middle");
                text.assign("dominant-baseline", "text-after-edge");
            }
        }

        let mut group = Group::new().add(circle).add(line).add(text);

        if let Some(class) = class {
            group.assign("class", class);
        } else {
            group.assign("class", "default");
        }
        if let Some(style) = style {
            group.assign("style", style);
        }
        self.annotations.append(group);
        self
    }

    /// Draw a Path using the Chart Coordinates onto the
    /// `scaled_layer``
    pub fn draw_path(mut self, mut path: Path, class: Option<&str>, style: Option<&str>) -> Self {
        if let Some(style) = style {
            path.assign("style", style);
        }
        if let Some(class) = class {
            path.assign("class", class);
        }
        self.xyplot.append(path);
        self
    }

    /// Add Data, composed of e.g. move_to and line_to operations, to the Chart
    /// using scaled coordinates.
    pub fn data(self, data: Data, class: Option<&str>, style: Option<&str>) -> Self {
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
    pub fn draw_shape(
        self,
        data: impl IntoIterator<Item = (f64, f64)>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        self.draw_path(to_path(data, true), class, style)
    }

    /*
    pub fn add_style(self, select: &str, style: &str) -> Self {
        let mut styles = self.svgdoc.styles.borrow_mut();
        styles.insert(select.to_string(), style.to_string());
        self
    }
     */
}

impl<'a> From<Chart<'a>> for Layer {
    fn from(chart: Chart) -> Self {
        let layer = Layer::new()
            .set("id", format!("xyplot-{}", last_id()))
            .set("class", "xyplot")
            .add(chart.xyplot)
            .add(chart.annotations)
            .add(chart.axes);
        // if !chart.canvas.styles.borrow().is_empty() {
        //     group = group.set("style", chart.canvas.styles.borrow().to_string());
        // }
        layer
    }
}

impl Display for Chart<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.xyplot.0)
    }
}

impl<'a> Deref for Chart<'a> {
    type Target = Group;

    fn deref(&self) -> &Self::Target {
        &self.xyplot.0
    }
}

impl<'a> DerefMut for Chart<'a> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.xyplot.0
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
