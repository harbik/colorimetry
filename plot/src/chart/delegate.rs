/// Macro to delegate all XYChart methods to a field of a struct, to avoid using Deref and DerefMut
/// Add methods here when adding new methods to XYChart
/// Usage: `delegate_xy_chart_methods!(MyStruct, chart_field);`
#[macro_export]
macro_rules! delegate_xy_chart_methods {
    ($struct_type:ty, $field:ident) => {
        impl $struct_type {
            // Basic drawing methods

            pub fn plot_grid(mut self, x_step: f64, y_step: f64, style_attr: StyleAttr) -> Self {
                self.$field = self.$field.plot_grid(x_step, y_step, style_attr);
                self
            }

            pub fn plot_poly_line(
                mut self,
                points: impl IntoIterator<Item = (f64, f64)>,
                style_attr: StyleAttr,
            ) -> Self {
                self.$field = self.$field.plot_poly_line(points, style_attr);
                self
            }

            pub fn plot_shape(
                mut self,
                points: impl IntoIterator<Item = (f64, f64)>,
                style_attr: StyleAttr,
            ) -> Self {
                self.$field = self.$field.plot_shape(points, style_attr);
                self
            }

            pub fn draw_data(
                mut self,
                layer: &str,
                data: svg::node::element::path::Data,
                style_attr: StyleAttr,
            ) -> Self {
                self.$field = self.$field.draw_data(layer, data, style_attr);
                self
            }

            pub fn draw_path(
                mut self,
                layer: &str,
                path: svg::node::element::Path,
                style_attr: StyleAttr,
            ) -> Self {
                self.$field = self.$field.draw_path(layer, path, style_attr);
                self
            }

            pub fn plot_image(
                mut self,
                image: impl Into<svg::node::element::Image>,
                style_attr: StyleAttr,
            ) -> Self {
                self.$field = self.$field.plot_image(image, style_attr);
                self
            }
            pub fn ticks(
                mut self,
                x_step: f64,
                y_step: f64,
                length: i32,
                style_attr: StyleAttr,
            ) -> Self {
                self.$field = self.$field.ticks(x_step, y_step, length, style_attr);
                self
            }

            pub fn x_labels(mut self, step: f64, offset: usize) -> Self {
                self.$field = self.$field.x_labels(step, offset);
                self
            }

            pub fn y_labels(mut self, step: f64, offset: usize) -> Self {
                self.$field = self.$field.y_labels(step, offset);
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

            // Annotation methods
            pub fn label_pin(
                mut self,
                cxy: (f64, f64),
                r: f64,
                angle_and_length: (i32, i32),
                text: impl AsRef<str>,
                style_attr: StyleAttr,
            ) -> Self {
                self.$field = self
                    .$field
                    .label_pin(cxy, r, angle_and_length, text, style_attr);
                self
            }
        }
    };
}
