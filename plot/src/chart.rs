mod delegate;

mod range;
pub use range::{ScaleRange, ScaleRangeIterator, ScaleRangeWithStep, ScaleValue};

use std::{collections::HashMap, ops::RangeBounds, rc::Rc};

use svg::{
    node::element::{path::Data, ClipPath, Definitions, Group, Image, Line, Path, Text, SVG},
    Node,
};

use crate::{
    layer::Layer, rendable::Rendable, round_to_default_precision, set_class_and_style,
    view::ViewParameters,
};

pub type CoordinateTransform = Rc<dyn Fn((f64, f64)) -> (f64, f64)>;

#[derive(Clone)]
pub struct XYChart {
    // inputs
    pub id: String,          // unique id for the chart
    pub x_range: ScaleRange, // min, max for x axis
    pub y_range: ScaleRange, // min, max for y axis
    pub plot_width: u32,     // width and height of the true plot area,
    pub plot_height: u32,    // excluding axis and margins

    // require inputs
    pub to_plot: CoordinateTransform,
    pub to_world: CoordinateTransform,

    // initially empty, filled by the chart methods
    pub view_parameters: ViewParameters,

    pub layers: HashMap<&'static str, Layer>,
    pub clip_paths: Vec<ClipPath>,
    pub margins: [i32; 4], // top, right, bottom, left
}

impl XYChart {
    pub const ANNOTATE_SEP: i32 = 2;
    pub const LABEL_HEIGHT: i32 = 16;
    pub const DESCRIPTION_HEIGHT: i32 = 20;
    pub const DESCRIPTION_SEP: i32 = 20;
    pub const DESCRIPTION_OFFSET: i32 = Self::LABEL_HEIGHT + Self::DESCRIPTION_SEP;

    pub fn new(
        id: impl AsRef<str>,
        plot_width_and_height: (u32, u32),
        ranges: (impl RangeBounds<f64>, impl RangeBounds<f64>),
        class_and_style: (Option<&str>, Option<&str>),
    ) -> XYChart {
        let id = id.as_ref().to_string();
        let (plot_width, plot_height) = plot_width_and_height;
        let (x_range, y_range) = (ScaleRange::new(ranges.0), ScaleRange::new(ranges.1));
        let (class, style) = class_and_style;

        // output coordinates with 0.1px precision
        let to_plot = Rc::new(move |xy: (f64, f64)| {
            let w = plot_width as f64;
            let h = plot_height as f64;
            world_to_plot_coordinates(xy.0, xy.1, &x_range, &y_range, w, h)
        });

        let to_world = Rc::new(move |xy: (f64, f64)| {
            let w = plot_width as f64;
            let h = plot_height as f64;
            plot_to_world_coordinates(xy.0, xy.1, &x_range, &y_range, w, h)
        });

        let mut clip_paths = Vec::new();

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
                .set("id", format!("clip-{id}"))
                .add(path.clone()),
        );
        let mut plot_layer = Layer::new();
        plot_layer.assign("clip-path", format!("url(#clip-{id})"));

        // don't add a background if no style or class is given
        if style.is_some() || class.is_some() {
            if let Some(style) = style {
                path = path.set("style", style);
            }
            if let Some(class) = class {
                path = path.set("class", class);
            }
            plot_layer.append(path);
        }

        // create the layers
        let mut layers = HashMap::new();
        layers.insert("plot", plot_layer);

        let mut annotations_layer = Layer::new();
        annotations_layer.assign("class", "annotations");
        layers.insert("annotations", annotations_layer);

        let mut axes_layer = Layer::new();
        axes_layer.assign("class", "axes");
        layers.insert("axes", axes_layer);

        let view_box = ViewParameters::new(0, 0, plot_width, plot_height, plot_width, plot_height);
        XYChart {
            id,
            view_parameters: view_box,
            plot_height,
            plot_width,
            x_range,
            y_range,
            layers,
            clip_paths,
            margins: [0i32; 4], // top, right, bottom, left
            to_plot,
            to_world,
        }
    }

    /// Adds ticks to all the sides of the plot.
    /// A negative length value produces inward ticks, and a positive value outward ticks.
    pub fn add_ticks(
        mut self,
        x_step: f64,
        y_step: f64,
        length: i32,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let mut data = Data::new();
        let to_plot = self.to_plot.clone();
        for x in self.x_range.iter_with_step(x_step) {
            let (px, py) = to_plot((x, self.y_range.start));
            data = data.move_to((px, py)).line_to((px, py + length as f64));

            let (px, py) = to_plot((x, self.y_range.end));
            data = data.move_to((px, py)).line_to((px, py - length as f64));
        }
        for y in self.y_range.iter_with_step(y_step) {
            let (px, py) = to_plot((self.x_range.start, y));
            data = data.move_to((px, py)).line_to((px - length as f64, py));

            let (px, py) = to_plot((self.x_range.end, y));
            data = data.move_to((px, py)).line_to((px + length as f64, py));
        }
        // Extend view coordinates if required
        if length > 0 {
            self.margins.iter_mut().for_each(|v| {
                if *v < length {
                    *v = length;
                }
            });
            self.update_view();
        }
        self.draw_data("axes", data, class, style)
    }
    pub fn add_x_labels(mut self, step: f64, offset: usize) -> Self {
        let range_with_step = ScaleRangeWithStep::new(self.x_range, step);
        let y = self.y_range.start;
        let to_plot = self.to_plot.clone();
        let mut x_labels = Group::new().set("class", "x-labels");
        for x in range_with_step.iter() {
            let display_value = format!("{}", ScaleValue(x, step));
            let (px, py) = to_plot((x, y));
            let txt = Text::new(display_value)
                .set("x", px)
                .set("y", py + offset as f64)
                .set("text-anchor", "middle")
                .set("dominant-baseline", "text-before-edge");
            x_labels.append(txt);
        }
        self.layers.get_mut("axes").unwrap().append(x_labels);
        self.margins[2] = self.margins[2].max(Self::LABEL_HEIGHT + offset as i32);
        self.update_view();
        self
    }
    pub fn add_y_labels(mut self, step: f64, offset: usize) -> Self {
        let range_with_step = ScaleRangeWithStep::new(self.y_range, step);
        let x = self.x_range.start;
        let to_plot = self.to_plot.clone();
        let mut y_labels = Group::new().set("class", "y-labels");
        for y in range_with_step.iter() {
            let display_value = format!("{}", ScaleValue(y, step));
            let (px, py) = to_plot((x, y));
            let txt = Text::new(display_value)
                .set("x", px - offset as f64)
                .set("y", py)
                .set("text-anchor", "middle")
                .set("dominant-baseline", "text-after-edge")
                .set(
                    "transform",
                    format!("rotate(-90, {}, {})", px - offset as f64, py),
                );
            y_labels.append(txt);
        }
        self.layers.get_mut("axes").unwrap().append(y_labels);
        self.margins[3] = self.margins[3].max(Self::LABEL_HEIGHT + offset as i32);
        self.update_view();
        self
    }

    pub fn x_axis_description(mut self, description: &str) -> Self {
        let x_middle = (self.x_range.start + self.x_range.end) / 2.0;
        let y = self.y_range.start;
        let (px, py) = (self.to_plot)((x_middle, y));
        let text = Text::new(description)
            .set("x", px)
            .set("y", py + Self::DESCRIPTION_OFFSET as f64)
            .set("text-anchor", "middle")
            .set("dominant-baseline", "text-before-edge")
            .set("class", "axis-description");
        self.layers.get_mut("axes").unwrap().append(text);
        self.margins[0] = self.margins[0].max(Self::DESCRIPTION_HEIGHT + Self::DESCRIPTION_OFFSET);
        self.update_view();
        self
    }

    /// Add a y-axis description, positioned to the left of the axis.
    ///
    /// The description is centered vertically along the y-axis.
    pub fn y_axis_description(mut self, description: &str) -> Self {
        let y_middle = (self.y_range.start + self.y_range.end) / 2.0;
        let x = self.x_range.start;
        let (px, py) = (self.to_plot)((x, y_middle));
        let text = Text::new(description)
            .set("x", px - Self::DESCRIPTION_OFFSET as f64)
            .set("y", py)
            .set("text-anchor", "middle")
            .set("dominant-baseline", "text-after-edge")
            .set("class", "axis-description")
            .set(
                "transform",
                format!(
                    "rotate(-90, {}, {})",
                    px - Self::DESCRIPTION_OFFSET as f64,
                    py
                ),
            );
        self.layers.get_mut("axes").unwrap().append(text);
        self.margins[3] = self.margins[3].max(Self::DESCRIPTION_HEIGHT + Self::DESCRIPTION_OFFSET);
        self.update_view();
        self
    }

    /*
    pub fn add_axis(
        mut self,
        description: Option<&str>,
        side: AxisSide,
        step: f64,
        tick_length: i32,
        show_labels: bool,
        class: Option<&str>,
        style: Option<&str>,
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
        let range_with_step = ScaleRangeWithStep::new(range, step);

        let x_axis = Axis::new(
            description,
            (0, 0, self.plot_width, self.plot_height),
            range_with_step,
            side,
            tick_length,
            show_labels,
            class,
        );

        self.layers.get_mut("axes").unwrap().append(Group::from(x_axis));
        self.update_view();
        self
    }
     */

    pub fn draw_grid(
        self,
        x_step: f64,
        y_step: f64,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let mut data = Data::new();
        let on_canvas = self.to_plot.clone();
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
        self.draw_data("plot", data, class, style)
    }

    pub fn draw_image(
        mut self,
        image: impl Into<Image>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let mut image: Image = image.into();
        image = set_class_and_style(image, class, style);
        self.layers.get_mut("plot").unwrap().append(image);
        self
    }

    pub fn annotate(
        mut self,
        cxy: (f64, f64),
        r: f64,
        angle_and_length: (i32, i32),
        text: impl AsRef<str>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let (angle, len) = angle_and_length;
        let angle = 360 - ((angle + 360) % 360);
        let (cx, cy) = cxy;
        let (cx, cy) = (self.to_plot)((cx, cy));
        let dx = len as f64 * (angle as f64).to_radians().cos();
        let dy = len as f64 * (angle as f64).to_radians().sin();
        let circle = svg::node::element::Circle::new()
            .set("cx", cx)
            .set("cy", cy)
            .set("r", r);
        let line = Line::new()
            .set("x1", cx)
            .set("y1", cy)
            .set("x2", cx + dx)
            .set("y2", cy + dy);

        let dxt = (len + Self::ANNOTATE_SEP) as f64 * (angle as f64).to_radians().cos();
        let dyt = (len + Self::ANNOTATE_SEP) as f64 * (angle as f64).to_radians().sin();
        let mut text = Text::new(text.as_ref())
            .set("x", cx + dxt)
            .set("y", cy + dyt);

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
        group = set_class_and_style(group, class, style);
        self.layers.get_mut("annotations").unwrap().append(group);
        self
    }

    /// Draw a Path using the Chart Coordinates onto the selected layer
    pub fn draw_path(
        mut self,
        layer: &str,
        mut path: Path,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        if let Some(layer) = self.layers.get_mut(layer) {
            path = set_class_and_style(path, class, style);
            layer.append(path);
        } else {
            panic!("unknown layer");
        }
        self
    }

    /// Add Data, composed of e.g. move_to and line_to operations, to chose layer
    pub fn draw_data(
        self,
        layer: &str,
        data: Data,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let path = Path::new().set("d", data);
        self.draw_path(layer, path, class, style)
    }

    /// Add a path to the chart, using coordinates from an iterator.
    pub fn draw_line(
        self,
        data: impl IntoIterator<Item = (f64, f64)>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let on_canvas = self.to_plot.clone();
        let iter_canvas = data.into_iter().map(|xy| (on_canvas)(xy));
        self.draw_path("plot", to_path(iter_canvas, false), class, style)
    }

    /// Add a path that connects the coordinates in an interator and closes the path.
    pub fn draw_shape(
        self,
        data: impl IntoIterator<Item = (f64, f64)>,
        class: Option<&str>,
        style: Option<&str>,
    ) -> Self {
        let on_canvas = self.to_plot.clone();
        let iter_canvas = data.into_iter().map(|xy| (on_canvas)(xy));
        self.draw_path("plot", to_path(iter_canvas, true), class, style)
    }

    pub fn update_view(&mut self) {
        let vx = -(self.margins[3]);
        let vy = -(self.margins[0]);
        // add margins[3] twice because of the left shift
        let vw = self.plot_width + self.margins[1] as u32 + 2 * self.margins[3] as u32;
        // add margins[0] twice because of the top shift
        let vh = self.plot_height + 2 * self.margins[0] as u32 + self.margins[2] as u32;

        self.view_parameters.set_view_box(vx, vy, vw, vh);
        self.set_width(vw);
        self.set_height(vh);
    }
}

/*
impl Display for XYChart {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.plot_layer.0)
    }
}
 */

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

/// Converts from world coordinates to plot coordinates.
/// The x and y coordinates are scaled to the plot width and height.
fn world_to_plot_coordinates(
    x: f64,
    y: f64,
    x_range: &ScaleRange,
    y_range: &ScaleRange,
    w: f64,
    h: f64,
) -> (f64, f64) {
    let x_canvas = x_range.scale(x) * w;
    let y_canvas = h - (y_range.scale(y) * h);
    (
        round_to_default_precision(x_canvas),
        round_to_default_precision(y_canvas),
    )
}

pub fn plot_to_world_coordinates(
    x: f64,
    y: f64,
    x_range: &ScaleRange,
    y_range: &ScaleRange,
    w: f64,
    h: f64,
) -> (f64, f64) {
    let x_world = x_range.unscale(x / w);
    let y_world = y_range.unscale_descent(y / h);
    (x_world, y_world)
}
#[test]
fn test_plot_to_world_coordinates() {
    use approx::assert_abs_diff_eq;
    let x_range = ScaleRange::new(0.0..=1.0);
    let y_range = ScaleRange::new(0.0..=1.0);
    let (x, y) = plot_to_world_coordinates(100.0, 200.0, &x_range, &y_range, 400.0, 300.0);
    assert_abs_diff_eq!(x, 0.25, epsilon = 1e-10);
    assert_abs_diff_eq!(y, 1.0 / 3.0, epsilon = 1e-10);
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
            .add(self.layers.get("axes").unwrap().clone())
            .add(self.layers.get("plot").unwrap().clone())
            .add(self.layers.get("annotations").unwrap().clone())
    }
}
