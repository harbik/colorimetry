use approx::AbsDiffEq;
use nalgebra::Vector3;

use crate::lab::CieLab;

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

    /// Sets the illuminance of the color represented by this `RelXYZ` instance.
    /// # Arguments
    /// - `illuminance`: The desired illuminance level in lux.
    ///
    /// # Returns
    /// A new `RelXYZ` instance with the XYZ values scaled to the specified illuminance.
    ///
    /// # Panics
    /// Panics if the Y value of the white point or the illuminance is less than or equal to zero.
    pub fn set_illuminance(mut self, illuminance: f64) -> Self {
        let yn = self.white_point.xyz.y;
        if yn > f64::EPSILON && illuminance > f64::EPSILON {
            let scale = illuminance / yn;
            self.xyz *= scale;
            self.white_point.xyz *= scale;
            self
        } else {
            panic!("Illuminance and Y of white point must be greater than zero");
        }
    }

    /// Checks if a related XYZ color is valid by converting it to CIELAB and back,
    /// and verifying the result is consistent, finite, and non-negative.
    ///
    /// # Returns
    /// - `true` if the XYZ value is perceptually valid and reversible.
    /// - `false` if the round-trip introduces a large error or the output contains invalid values.
    ///
    /// # Why This Works
    /// The XYZ → CIELAB → XYZ transformation is only reliable for physically meaningful colors.
    /// If the input XYZ is too far from the reference white, or contains negative components,
    /// the LAB model may produce invalid results or large reversibility errors.
    pub fn is_valid(&self) -> bool {
        let lab = CieLab::from_xyz(*self);
        let xyz_back = lab.xyz();
        let same = self.abs_diff_eq(&xyz_back, 1E-5);
        same && xyz_back.values()[0]
            .into_iter()
            .all(|v| v >= 0.0 && v.is_finite())
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observer::Observer::Cie1931;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_new() {
        let xyz = [1.0, 2.0, 3.0];
        let white = XYZ::new([100.0, 100.0, 100.0], Cie1931);
        let rel_xyz = RelXYZ::new(xyz, white);
        assert_abs_diff_eq!(rel_xyz.xyz.x, 1.0);
        assert_abs_diff_eq!(rel_xyz.xyz.y, 2.0);
        assert_abs_diff_eq!(rel_xyz.xyz.z, 3.0);
        assert_abs_diff_eq!(rel_xyz.white_point.xyz.x, 100.0);
    }

    #[test]
    fn test_from_vec() {
        let xyz = Vector3::new(1.0, 2.0, 3.0);
        let white = XYZ::new([100.0, 100.0, 100.0], Cie1931);
        let rel_xyz = RelXYZ::from_vec(xyz, white);
        assert_abs_diff_eq!(rel_xyz.xyz.x, 1.0);
        assert_abs_diff_eq!(rel_xyz.xyz.y, 2.0);
        assert_abs_diff_eq!(rel_xyz.xyz.z, 3.0);
    }

    #[test]
    fn test_with_d65() {
        let xyz = XYZ::new([1.0, 2.0, 3.0], Cie1931);
        let rel_xyz = RelXYZ::with_d65(xyz);
        assert_abs_diff_eq!(rel_xyz.xyz.x, 1.0);
        assert_eq!(rel_xyz.white_point, Cie1931.xyz_d65());
    }

    #[test]
    fn test_with_d50() {
        let xyz = XYZ::new([1.0, 2.0, 3.0], Cie1931);
        let rel_xyz = RelXYZ::with_d50(xyz);
        assert_abs_diff_eq!(rel_xyz.xyz.x, 1.0);
        assert_eq!(rel_xyz.white_point, Cie1931.xyz_d50());
    }

    #[test]
    fn test_set_illuminance() {
        let white = XYZ::new([50.0, 50.0, 50.0], Cie1931);
        let rel_xyz = RelXYZ::new([5.0; 3], white).set_illuminance(100.0);
        assert_abs_diff_eq!(rel_xyz.xyz.y, 10.0);
        assert_abs_diff_eq!(rel_xyz.white_point.xyz.y, 100.0);
    }

    #[test]
    #[should_panic]
    fn test_set_illuminance_zero() {
        let xyz = XYZ::new([1.0, 2.0, 3.0], Cie1931);
        let white = XYZ::new([100.0, 100.0, 100.0], Cie1931);
        let rel_xyz = RelXYZ::from_xyz(xyz, white).unwrap();
        rel_xyz.set_illuminance(0.0);
        // This should panic because illuminance is zero
    }

    #[test]
    fn test_spectral_locus_round_trip() {
        use crate::{illuminant::CieIlluminant, observer::Observer::Cie1931};

        let sl = Cie1931.spectral_locus(CieIlluminant::D65);
        for (_w, rxyz) in sl {
            let lab = CieLab::from_xyz(rxyz);
            let xyz_back = lab.xyz();
            approx::assert_abs_diff_eq!(rxyz, xyz_back, epsilon = 1E-6)
        }
    }

    #[test]
    fn test_spectral_locus_round_trip_print() {
        use crate::{illuminant::CieIlluminant, observer::Observer::Cie1931};

        let sl = Cie1931.spectral_locus(CieIlluminant::D65);
        for (w, rxyz) in sl {
            print!("{}, {:.4?}", w, rxyz.values()[0]);
            let lab = CieLab::from_xyz(rxyz);
            let xyz_back = lab.xyz();
            println!("{:.4?}", xyz_back.values()[0]);
        }
    }
}
