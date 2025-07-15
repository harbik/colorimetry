use nalgebra::{Matrix3, Vector3};

use crate::axis::ChartRange;

// SVG transform matrix: matrix(a, b, c, d, e, f)
// which is: [a c e]
//           [b d f]
//           [0 0 1]
#[derive(Debug, Clone)]
/// A struct to handle transformations between chart coordinates and canvas coordinates.
pub struct CoordinateTransform {
    to_canvas_matrix: Matrix3<f64>,
    to_scaled_matrix: Matrix3<f64>,
}

impl CoordinateTransform {
    pub fn new(target: [u32; 4], range_x: ChartRange, range_y: ChartRange) -> Self {
        let [left, top, width, height] = target;
        let scale_x = width as f64 / range_x.span();
        let scale_y = height as f64 / range_y.span();

        // SVG y-axis increases downward, so invert y scaling
        let scale_y = -scale_y;
        let translate_x = left as f64 - range_x.start * scale_x;
        let translate_y = top as f64 + height as f64 + range_y.start * scale_y;

        let to_canvas_matrix = Matrix3::new(
            scale_x,
            0.0,
            translate_x,
            0.0,
            scale_y,
            translate_y,
            0.0,
            0.0,
            1.0,
        );

        let to_scale_matrix = to_canvas_matrix.try_inverse().expect("Matrix is not invertible");

        CoordinateTransform {
            to_canvas_matrix,
            to_scaled_matrix: to_scale_matrix,
        }
    }


    pub fn on_canvas(
        &self,
        x: f64,
        y: f64,
    ) -> (f64, f64) {
        let point = Vector3::new(x, y, 1.0);
        let transformed = self.to_canvas_matrix * point;
        (transformed[(0, 0)], transformed[(1, 0)])
    }

    pub fn scaled(
        &self,
        h: u32,
        v: u32,
    ) -> (f64, f64) {
        let point = Vector3::new(h as f64, v as f64, 1.0);
        let &[x,y, _] : &[f64;3] = (self.to_scaled_matrix * point).as_ref();
        (x, y)

    }

 
    pub fn to_chart_string(&self) -> String {
        format!(
            "matrix({} {} {} {} {} {})",
            self.to_canvas_matrix[(0, 0)],
            self.to_canvas_matrix[(1, 0)],
            self.to_canvas_matrix[(0, 1)],
            self.to_canvas_matrix[(1, 1)],
            self.to_canvas_matrix[(0, 2)],
            self.to_canvas_matrix[(1, 2)]
        )
    }

    pub fn to_canvas_string(&self) -> String {
        format!(
            "matrix({} {} {} {} {} {})",
            self.to_scaled_matrix[(0, 0)],
            self.to_scaled_matrix[(1, 0)],
            self.to_scaled_matrix[(0, 1)],
            self.to_scaled_matrix[(1, 1)],
            self.to_scaled_matrix[(0, 2)],
            self.to_scaled_matrix[(1, 2)]
        )
    }
}

