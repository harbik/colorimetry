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
                style: Option<&str>,
            ) -> Self {
                self.$field = self.$field.add_ticks(x_step, y_step, length, class, style);
                self
            }

            pub fn add_x_labels(mut self, step: f64, offset: usize) -> Self {
                self.$field = self.$field.add_x_labels(step, offset);
                self
            }

            pub fn add_y_labels(mut self, step: f64, offset: usize) -> Self {
                self.$field = self.$field.add_y_labels(step, offset);
                self
            }

            pub fn x_axis_description(mut self, description: &str) -> Self {
                self.$field = self.$field.x_axis_description(description);
                self
            }

            pub fn y_axis_description(mut self, description: &str) -> Self {
                self.$field = self.$field.y_axis_description(description);
                self
            }

            /*
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
             */

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
