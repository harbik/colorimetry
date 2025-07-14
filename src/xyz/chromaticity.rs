use nalgebra::Vector2;

/// A chromaticity coordinate with x and y values.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Chromaticity {
    xy: Vector2<f64>,
}

impl Chromaticity {
    /// Returns a new `Chromaticity` object with the given x and y coordinates.
    pub const fn new(x: f64, y: f64) -> Self {
        let xy = Vector2::new(x, y);
        Self { xy }
    }

    /// Returns the x coordinate.
    pub fn x(self) -> f64 {
        self.xy.x
    }

    /// Returns the y coordinate.
    pub fn y(self) -> f64 {
        self.xy.y
    }

    /// Returns the chromaticity coordinates as an array in the format [x, y].
    ///
    /// ```
    /// # use colorimetry::xyz::Chromaticity;
    /// let [x, y] = [0.3127, 0.3290];
    /// let chromaticity = Chromaticity::new(x, y);
    /// assert_eq!(chromaticity.to_array(), [x, y]);
    /// ```
    pub fn to_array(self) -> [f64; 2] {
        *self.xy.as_ref()
    }

    /// Returns the chromaticity coordinate as an `nalgebra` vector.
    pub const fn to_vector(self) -> Vector2<f64> {
        self.xy
    }

    pub fn to_tuple(&self) -> (f64, f64) {
        (self.xy.x, self.xy.y)
    }
}

impl Default for Chromaticity {
    /// Returns the default chromaticity coordinate, which is the Equal Energy  white point
    fn default() -> Self {
        Self::new(1.0 / 3.0, 1.0 / 3.0)
    }
}
