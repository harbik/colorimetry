use std::{
    fmt::Display,
    ops::{Deref, DerefMut, RangeBounds},
    rc::Rc,
};

use svg::{
    node::element::{path::Data, ClipPath, Definitions, Group, Line, Path, Text, SVG},
    Node,
};

use crate::{
    assign_class_and_style, axis::{Axis, AxisSide, ChartRange}, layer::Layer, round_to_default_precision, rendable::Rendable,  view::ViewParameters
};

#[derive(Clone)]
pub struct XYChart {
    pub id: String,          // unique id for the chart
    pub view_parameters: ViewParameters,
    pub x_range: ChartRange, // min, max for x axis
    pub y_range: ChartRange, // min, max for y axis
    pub plot: Layer,
    pub axes: Layer,
    pub annotations: Layer,
    pub clip_paths: Vec<ClipPath>,
    pub margins: [i32; 4],   // top, right, bottom, left
    pub on_canvas: Rc<dyn Fn((f64, f64)) -> (f64, f64)>,
    pub plot_width: u32,
    pub plot_height: u32,
}

impl XYChart {
    pub const ANNOTATE_SEP: i32 = 2;
    pub const LABEL_HEIGHT: i32 = 16;
    pub const DESCRIPTION_HEIGHT: i32 = 20;
    pub const DESCRIPTION_SEP: i32 = 15;
    pub const DESCRIPTION_OFFSET: i32 = Self::LABEL_HEIGHT + Self::DESCRIPTION_SEP;

    pub fn new(
        id: impl AsRef<str>,
        plot_width: u32,
        plot_height: u32,
        x_range: impl RangeBounds<f64>,
        y_range: impl RangeBounds<f64>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> XYChart {
        let id = id.as_ref().to_string();
        let x_range = ChartRange::new(x_range);
        let y_range = ChartRange::new(y_range);

        // output coordinates with 0.1px precision
        let on_canvas = Rc::new(move |xy: (f64, f64)| {
            let w = plot_width as f64;
            let h = plot_height as f64;
            convert_to_plot_coordinates(xy.0, xy.1, &x_range, &y_range, w, h)
        });
        let mut clip_paths = Vec::new();
        let mut plot = Layer::new();
        let mut path = to_path(
            [
                (0f64, 0f64),
                (plot_width as f64, 0f64),
                (plot_width as f64, plot_height as f64),
                (0f64, plot_height as f64),
            ],
            true,
        );
        // chart always gets a clip path
        clip_paths.push(
            ClipPath::new()
                .set("id", format!("clip-{}", id))
                .add(path.clone()),
        );
        plot = plot.set("clip-path", format!("url(#clip-{})", id));
        // don't add a background if no style or class is given
        if style.is_some() || class.is_some() {
            if let Some(style) = style {
                path = path.set("style", style);
            }
            if let Some(class) = class {
                path = path.set("class", class);
            }
            plot.append(path);
        }
        let annotations = Layer::new()
            .set("id", format!("annotations-{}", id))
            .set("class", "annotations");

        let axes = Layer::new()
            .set("id", format!("axes-{}", id))
            .set("class", "axes");
        let view_box = ViewParameters::new(0, 0, plot_width, plot_height, plot_width, plot_height);
        XYChart {
            id,
            view_parameters: view_box,
            plot_height,
            plot_width,
            x_range,
            y_range,
            plot,
            annotations,
            axes,
            clip_paths,
            margins: [0i32; 4], // top, right, bottom, left
            on_canvas,
            
        }
    }

    pub fn add_axis(
        mut self,
        description: Option<&str>,
        side: AxisSide,
        step: f64,
        tick_length: i32,
        show_labels: bool,
        class: Option<&str>,
    ) -> Self {
        let mut margin = if show_labels {
            tick_length + Self::LABEL_HEIGHT
        } else {
            tick_length
        };

        if description.is_some() {
            margin += Self::DESCRIPTION_HEIGHT + Self::DESCRIPTION_OFFSET;
        };

        match side {
            AxisSide::Top => {
                self.margins[0] = self.margins[0].max(margin); // top margin
            }
            AxisSide::Right => {
                self.margins[1] = self.margins[1].max(margin); // right margin
            }
            AxisSide::Bottom => {
                self.margins[2] = self.margins[2].max(margin); // bottom margin
            }
            AxisSide::Left => {
                self.margins[3] = self.margins[3].max(margin); // left margin
            }
        }
        let range = match side {
            AxisSide::Bottom | AxisSide::Top => self.x_range,
            AxisSide::Left | AxisSide::Right => self.y_range,
        };
        let x_axis = Axis::new(
            description,
            (0, 0, self.view_parameters.width(), self.view_parameters.height()),
            range,
            step,
            side,
            tick_length,
            show_labels,
            class,
        );
        self.axes.append(Group::from(x_axis));
        self.update_view();
        self
    }

    pub fn draw_grid(
        self,
        x_step: f64,
        y_step: f64,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let mut data = Data::new();
        let on_canvas = self.on_canvas.clone();
        for x in self.x_range.iter_with_step(x_step) {
            data = data
                .move_to(on_canvas((x, self.y_range.start)))
                .line_to(on_canvas((x, self.y_range.end)));
        }
        for y in self.y_range.iter_with_step(y_step) {
            data = data
                .move_to(on_canvas((self.x_range.start, y)))
                .line_to(on_canvas((self.x_range.end, y)));
        }
        self.draw_data(data, class, style)
    }

    pub fn annotate(
        mut self,
        cxy: (f64, f64),
        r: f64,
        len: i32,
        angle: i32,
        text: impl AsRef<str>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let angle = 360 - ((angle + 360) % 360) as i32;
        let (cx, cy) = cxy;
        let (cx, cy) = (self.on_canvas)((cx, cy));
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
            0..55 | 305..=360 => {
                // easat
                text.assign("text-anchor", "start");
                text.assign("dominant-baseline", "middle");
            }
            55..125 => {
                // north
                text.assign("text-anchor", "middle");
                text.assign("dominant-baseline", "text-before-edge");
            }
            125..235 => {
                // west
                text.assign("text-anchor", "end");
                text.assign("dominant-baseline", "middle");
            }
            _ => {
                // south
                text.assign("text-anchor", "middle");
                text.assign("dominant-baseline", "text-after-edge");
            }
        }

        let mut group = Group::new().add(circle).add(line).add(text);
        assign_class_and_style(&mut group, class, style);
        self.annotations.append(group);
        self
    }

    /// Draw a Path using the Chart Coordinates onto the
    /// `scaled_layer``
    pub fn draw_path(mut self, mut path: Path, class: Option<&str>, style: Option<&str>) -> Self {
        assign_class_and_style(&mut path, class, style);
        self.plot.append(path);
        self
    }

    /// Add Data, composed of e.g. move_to and line_to operations, to the Chart
    /// using scaled coordinates.
    pub fn draw_data(self, data: Data, class: Option<&str>, style: Option<&str>) -> Self {
        let path = Path::new().set("d", data);
        self.draw_path(path, class, style)
    }

    /// Add a path to the chart, using coordinates from an iterator.
    pub fn draw_line(
        self,
        data: impl IntoIterator<Item = (f64, f64)>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let on_canvas = self.on_canvas.clone();
        let iter_canvas = data.into_iter().map(|xy| (on_canvas)(xy));
        self.draw_path(to_path(iter_canvas, false), class, style)
    }

    /// Add a path that connects the coordinates in an interator and closes the path.
    pub fn draw_shape(
        self,
        data: impl IntoIterator<Item = (f64, f64)>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let on_canvas = self.on_canvas.clone();
        let iter_canvas = data.into_iter().map(|xy| (on_canvas)(xy));
        self.draw_path(to_path(iter_canvas, true), class, style)
    }

    pub fn update_view(&mut self) {
        let vx = -(self.margins[3] as i32);
        let vy = -(self.margins[0] as i32);
        // add margins[3] twice because of the left shift
        let vw = self.plot_width + self.margins[1] as u32 + 2 * self.margins[3] as u32;
        // add margins[0] twice because of the top shift
        let vh = self.plot_height + 2 * self.margins[0] as u32 + self.margins[2] as u32;

        self.view_parameters.set_view_box(vx, vy, vw, vh);
        self.set_width(vw);
        self.set_height(vh);

    }


}

impl Display for XYChart {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.plot.0)
    }
}

impl Deref for XYChart {
    type Target = Group;

    fn deref(&self) -> &Self::Target {
        &self.plot.0
    }
}

impl DerefMut for XYChart {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.plot.0
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
}

fn convert_to_plot_coordinates(x: f64, y: f64, x_range: &ChartRange, y_range: &ChartRange, w: f64, h: f64) -> (f64, f64) {
    let x_canvas = x_range.scale(x) * w;
    let y_canvas = h - (y_range.scale(y) * h);
    (
        round_to_default_precision(x_canvas),
        round_to_default_precision(y_canvas),
    )
}

impl From<XYChart> for SVG {
    fn from(chart: XYChart) -> Self {
        chart.render()
    }
}

impl Rendable for XYChart {
    // required parameters
    fn view_parameters(&self) -> ViewParameters {
        self.view_parameters.clone()
    }

    fn set_view_parameters(&mut self, view_box: ViewParameters) {
        self.view_parameters = view_box;
    }

    fn render(&self) -> SVG {
        println!("Converting XYChart to SVG with id: {}", self.id);
        let mut defs = Definitions::new();
        for clip in self.clip_paths.iter() {
            defs.append(clip.clone());
        }
        SVG::new()
            .set("id", self.id.clone())
            .set("width", self.width())
            .set("height", self.height())
            .set("viewBox", self.view_parameters().view_box_str())
            .add(defs)
            .add(self.axes.clone())
            .add(self.plot.clone())
            .add(self.annotations.clone())
    }
    
}