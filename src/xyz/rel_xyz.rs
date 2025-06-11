use approx::AbsDiffEq;
use nalgebra::Vector3;

use super::XYZ;

/// # Related Tristimulus Values
///
/// Tristimulus Values for a given sample and reference white,
/// used to represent related colors as used in various color
/// models. Typically the reference white is normalized to have
/// an Y-value of 100
#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct RelXYZ {
    xyz: Vector3<f64>,
    white_point: XYZ,
}

impl RelXYZ {
    /// Creates a new `RelXYZ` instance with the given XYZ values and white point.
    /// # Arguments
    /// - `xyz`: A 3-element array representing the XYZ tristimulus values.
    /// - `white_point`: The reference white point as an `XYZ` instance.
    ///
    /// # Returns
    /// A new `RelXYZ` instance initialized with the provided XYZ values and white point.   
    pub fn new(xyz: [f64; 3], white_point: XYZ) -> Self {
        RelXYZ {
            xyz: xyz.into(),
            white_point,
        }
    }

    /// Creates a new `RelXYZ` instance with the given XYZ values and white point.
    /// # Arguments
    /// - `xyz`: A 3-element vector representing the XYZ tristimulus values.
    /// - `white_point`: The reference white point as an `XYZ` instance.
    ///
    /// # Returns
    /// A new `RelXYZ` instance initialized with the provided XYZ values and white point.   
    pub fn from_vec(xyz: Vector3<f64>, white_point: XYZ) -> Self {
        RelXYZ { xyz, white_point }
    }

    /// Creates a new `RelXYZ` instance from an `XYZ` instance and a reference white point.
    ///
    /// - `xyz`: An `XYZ` instance representing the color to be transformed.
    /// - `white_point`: An `XYZ` instance representing the reference white point.
    ///
    /// # Returns
    /// A new `RelXYZ` instance initialized with the XYZ values from the provided `XYZ` instance and the reference white point.
    /// # Errors
    /// Returns an error if the observer of the `xyz` and `white_point` do not match.   
    pub fn from_xyz(xyz: XYZ, white_point: XYZ) -> Result<Self, crate::Error> {
        if xyz.observer != white_point.observer {
            Err(crate::Error::RequireSameObserver)
        } else {
            Ok(RelXYZ {
                xyz: xyz.xyz,
                white_point,
            })
        }
    }

    /// Creates a new `RelXYZ` instance with the given XYZ values and a D65 white point.
    /// # Arguments
    /// - `xyz`: a XYZ tristimulus value.
    ///
    /// # Returns
    /// A new `RelXYZ` instance initialized with the provided XYZ value and a D65 white point.
    pub fn with_d65(xyz: XYZ) -> Self {
        let white_point = xyz.observer.xyz_d65();
        RelXYZ {
            xyz: xyz.xyz,
            white_point,
        }
    }

    /// Creates a new `RelXYZ` instance with the given XYZ values and a D50 white point.
    /// # Arguments
    /// - `xyz`: a XYZ tristimulus value.
    ///
    /// # Returns
    /// A new `RelXYZ` instance initialized with the provided XYZ value and a D50 white point.
    pub fn with_d50(xyz: XYZ) -> Self {
        let white_point = xyz.observer.xyz_d50();
        RelXYZ {
            xyz: xyz.xyz,
            white_point,
        }
    }

    /// Returns the XYZ tristimulus values of the color represented by this `RelXYZ` instance.
    pub fn xyz(&self) -> XYZ {
        XYZ::from_vecs(self.xyz, self.white_point.observer)
    }

    /// Returns the reference white point of this `RelXYZ` instance.
    pub fn white_point(&self) -> XYZ {
        self.white_point
    }

    /// Returns the XYZ tristimulus values of the color and the reference white point as a 2D array.
    ///
    /// The first row contains the XYZ values of the color, and the second row contains the XYZ values of the reference white point.        
    pub fn values(&self) -> [[f64; 3]; 2] {
        [self.xyz.into(), self.white_point.xyz.into()]
    }
}

impl AbsDiffEq for RelXYZ {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let xyz_eq = self.xyz.abs_diff_eq(&other.xyz, epsilon);
        let xyzn_eq = self.white_point.abs_diff_eq(&other.white_point, epsilon);
        xyz_eq && xyzn_eq
    }
}
