use nalgebra::Matrix3;

// SVG transform matrix: matrix(a, b, c, d, e, f)
// which is: [a c e]
//           [b d f]
//           [0 0 1]
pub struct TransformMatrix {
    to_chart: Matrix3<f64>,
    to_canvas: Matrix3<f64>,
}


impl TransformMatrix {
    pub fn new(target: [u32; 4], scale: [[f64; 2]; 2]) -> Self {
        let [[x_min, x_max], [y_min, y_max]] = scale;
        let [left, top, width, height] = target;
        let scale_x = width as f64 / (x_max - x_min);
        let scale_y = height as f64 / (y_max - y_min);

        // SVG y-axis increases downward, so invert y scaling
        let scale_y = -scale_y;
        let translate_x = left as f64 - x_min * scale_x;
        let translate_y = top as f64 + height as f64 + y_min * scale_y;

        let to_chart = Matrix3::new(
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

        let to_canvas = to_chart.try_inverse().expect("Matrix is not invertible");

        TransformMatrix {
            to_chart,
            to_canvas,
        }
    }

    pub fn to_chart_string(&self) -> String {
        format!(
            "matrix({} {} {} {} {} {})",
            self.to_chart[(0, 0)],
            self.to_chart[(1, 0)],
            self.to_chart[(0, 1)],
            self.to_chart[(1, 1)],
            self.to_chart[(0, 2)],
            self.to_chart[(1, 2)]
        )
    }

    pub fn to_canvas_string(&self) -> String {
        format!(
            "matrix({} {} {} {} {} {})",
            self.to_canvas[(0, 0)],
            self.to_canvas[(1, 0)],
            self.to_canvas[(0, 1)],
            self.to_canvas[(1, 1)],
            self.to_canvas[(0, 2)],
            self.to_canvas[(1, 2)]
        )
    }
}
