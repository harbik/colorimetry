
use std::{fmt::Display, ops::RangeBounds, rc::Rc};

use svg::{
    node::element::{path::Data, ClipPath, Definitions, Group, Image, Line, Path, Text, SVG},
    Node,
};

use crate::{
    axis::{Axis, AxisSide, ScaleRange, ScaleRangeWithStep},
    layer::Layer,
    rendable::Rendable,
    round_to_default_precision, set_class_and_style,
    view::ViewParameters,
};

pub type CoordinateTransform = Rc<dyn Fn((f64, f64)) -> (f64, f64)>;

#[derive(Clone)]
pub struct XYChart {
    // inputs
    pub id: String,          // unique id for the chart
    pub x_range: ScaleRange, // min, max for x axis
    pub y_range: ScaleRange, // min, max for y axis
    pub plot_width: u32,
    pub plot_height: u32,

    // require inputs
    pub to_plot: CoordinateTransform,
    pub to_world: CoordinateTransform,

    // initially empty, filled by the chart methods
    pub view_parameters: ViewParameters,
    
    // todo: make layers a hashmap: <"String", "Layer">
    pub plot_layer: Layer,
    pub axes_layer: Layer,
    pub annotations_layer: Layer,
    pub clip_paths: Vec<ClipPath>,
    pub margins: [i32; 4], // top, right, bottom, left
}

impl XYChart {
    pub const ANNOTATE_SEP: i32 = 2;
    pub const LABEL_HEIGHT: i32 = 16;
    pub const DESCRIPTION_HEIGHT: i32 = 20;
    pub const DESCRIPTION_SEP: i32 = 15;
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
        let mut plot_layer = Layer::new();

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
        plot_layer = plot_layer.set("clip-path", format!("url(#clip-{id})"));

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

        let annotations_layer = Layer::new()
            .set("id", format!("annotations-{id}"))
            .set("class", "annotations");

        let axes_layer = Layer::new()
            .set("id", format!("axes-{id}"))
            .set("class", "axes");
        let view_box = ViewParameters::new(0, 0, plot_width, plot_height, plot_width, plot_height);
        XYChart {
            id,
            view_parameters: view_box,
            plot_height,
            plot_width,
            x_range,
            y_range,
            plot_layer,
            annotations_layer,
            axes_layer,
            clip_paths,
            margins: [0i32; 4], // top, right, bottom, left
            to_plot,
            to_world,
        }
    }

    /// Adds ticks to all the sides of the plot.
    /// A negative length value produces inward ticks, and a positive value outward ticks.
    pub fn add_ticks(mut self, x_step: f64, y_step: f64, length: i32, class: Option<&str>, style: Option<&str>) -> Self {
        let mut data = Data::new();
        let to_plot = self.to_plot.clone();
        for x in self.x_range.iter_with_step(x_step) {

            let (px, py) = to_plot((x, self.y_range.start));
            data = data
                .move_to((px, py))
                .line_to((px, py + length as f64));

            let (px, py) = to_plot((x, self.y_range.end));
            data = data
                .move_to((px, py))
                .line_to((px, py - length as f64));
        }
        for y in self.y_range.iter_with_step(y_step) {

            let (px, py) = to_plot((self.x_range.start, y));
            data = data
                .move_to((px, py))
                .line_to((px - length as f64, py));

            let (px, py) = to_plot((self.x_range.end, y));
            data = data
                .move_to((px, py))
                .line_to((px + length as f64, py));
        }
        // Extend view coordinates if required
        if length > 0 {
            self.margins.iter_mut().for_each(|v|{
                if *v < length {
                    *v = length;
                }
            });
            self.update_view();
        }
        self.draw_data("axes", data, class, style)
        
    }

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
        
        self.axes_layer.append(Group::from(x_axis));
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
        self.plot_layer.append(image);
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
        self.annotations_layer.append(group);
        self
    }

    /// Draw a Path using the Chart Coordinates onto the selected layer
    pub fn draw_path(mut self, layer: &str,  mut path: Path, class: Option<&str>, style: Option<&str>) -> Self {
        path = set_class_and_style(path, class, style);
        match layer {
            "axes" => self.axes_layer.append(path),
            "plot" => self.plot_layer.append(path),
            "annotations" => self.annotations_layer.append(path),
            _ => panic!("layer must be one of 'axes' 'plot', or 'annotations'")
        }
        self
    }

    /// Add Data, composed of e.g. move_to and line_to operations, to chose layer
    pub fn draw_data(self, layer: &str, data: Data, class: Option<&str>, style: Option<&str>) -> Self {
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

impl Display for XYChart {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.plot_layer.0)
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
            .add(self.axes_layer.clone())
            .add(self.plot_layer.clone())
            .add(self.annotations_layer.clone())
    }
}

/// Macro to delegate all XYChart methods to a field of a struct, to avoid using Deref and DerefMut
/// Add methods here when adding new methods to XYChart
/// Usage: `delegate_xy_chart_methods!(MyStruct, chart_field);`
#[macro_export]
macro_rules! delegate_xy_chart_methods {
    ($struct_type:ty, $field:ident) => {
        impl $struct_type {
            // Basic drawing methods
            pub fn draw_grid(
                mut self,
                x_step: f64,
                y_step: f64,
                class: Option<&str>,
                style: Option<&str>,
            ) -> Self {
                self.$field = self.$field.draw_grid(x_step, y_step, class, style);
                self
            }

            pub fn draw_line(
                mut self,
                points: impl IntoIterator<Item = (f64, f64)>,
                class: Option<&str>,
                style: Option<&str>,
            ) -> Self {
                self.$field = self.$field.draw_line(points, class, style);
                self
            }

            pub fn draw_shape(
                mut self,
                points: impl IntoIterator<Item = (f64, f64)>,
                class: Option<&str>,
                style: Option<&str>,
            ) -> Self {
                self.$field = self.$field.draw_shape(points, class, style);
                self
            }

            pub fn draw_data(
                mut self,
                layer: &str,
                data: svg::node::element::path::Data,
                class: Option<&str>,
                style: Option<&str>,
            ) -> Self {
                self.$field = self.$field.draw_data(layer, data, class, style);
                self
            }

            pub fn draw_path(
                mut self,
                layer: &str,
                path: svg::node::element::Path,
                class: Option<&str>,
                style: Option<&str>,
            ) -> Self {
                self.$field = self.$field.draw_path(layer, path, class, style);
                self
            }

            pub fn draw_image(
                mut self,
                image: impl Into<svg::node::element::Image>,
                class: Option<&str>,
                style: Option<&str>,
            ) -> Self {
                self.$field = self.$field.draw_image(image, class, style);
                self
            }
            pub fn add_ticks(
                mut self, 
                x_step: f64, 
                y_step: f64, 
                length: i32, 
                class: Option<&str>, 
                style: Option<&str>
            ) -> Self {
                self.$field =
                    self.$field
                        .add_ticks(x_step, y_step, length, class, style);
                self

            }
            // Axis methods
            pub fn add_axis(
                mut self,
                description: Option<&str>,
                side: $crate::axis::AxisSide,
                step: f64,
                tick_length: i32,
                show_labels: bool,
                class: Option<&str>,
                style: Option<&str>
            ) -> Self {
                self.$field =
                    self.$field
                        .add_axis(description, side, step, tick_length, show_labels, class, style);
                self
            }

            // Annotation methods
            pub fn annotate(
                mut self,
                cxy: (f64, f64),
                r: f64,
                angle_and_length: (i32, i32),
                text: impl AsRef<str>,
                class: Option<&str>,
                style: Option<&str>,
            ) -> Self {
                self.$field = self
                    .$field
                    .annotate(cxy, r, angle_and_length, text, class, style);
                self
            }
        }
    };
}
