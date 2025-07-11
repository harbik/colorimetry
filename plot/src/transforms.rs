use nalgebra::{Matrix3, Vector3};

// SVG transform matrix: matrix(a, b, c, d, e, f)
// which is: [a c e]
//           [b d f]
//           [0 0 1]
#[derive(Debug, Clone)]
/// A struct to handle transformations between chart coordinates and canvas coordinates.
pub struct TransformMatrix {
    to_canvas_matrix: Matrix3<f64>,
    to_scaled_matrix: Matrix3<f64>,
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

        TransformMatrix {
            to_canvas_matrix,
            to_scaled_matrix: to_scale_matrix,
        }
    }

    pub fn canvas(
        &self,
        x: f64,
        y: f64,
    ) -> (u32, u32) {
        let point = Vector3::new(x, y, 1.0);
        let transformed = self.to_canvas_matrix * point;
        (transformed[(0, 0)] as u32, transformed[(1, 0)] as u32)
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


#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_transform_matrix() {
        let target = [0, 0, 100, 100];
        let scale = [[0.0, 1.0], [0.0, 1.0]];
        let matrix = TransformMatrix::new(target, scale);
        
        assert_eq!(matrix.to_chart_string(), "matrix(100 0 0 -100 0 100)");
        assert_eq!(matrix.to_canvas_string(), "matrix(0.01 0 -0 -0.01 -0 1)");
        
        let (h, v) = matrix.canvas(0.5, 0.5);
        assert_eq!((h, v), (50, 50));

        let (h, v) = matrix.canvas(0.0, 0.0);
        assert_eq!((h, v), (0, 100));
        
        let (h, v) = matrix.canvas(1.0, 1.0);
        assert_eq!((h, v), (100, 0));

        let (x, y) = matrix.scaled(50, 50);
        assert_eq!((x, y), (0.5, 0.5));

        let (x, y) = matrix.scaled(100, 0);
        assert_eq!((x, y), (1.0, 1.0));
    }
}