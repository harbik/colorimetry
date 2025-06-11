//! # CIECAM02 Color Appearance Model
//!
//! This module implements the **CIECAM02** appearance model.
//! It converts CIE XYZ tristimulus values into perceptual correlates of lightness (J),  
//! chroma (C), and hue angle (h) under specified viewing conditions, and supports:
//! - **Raw “Jab”** coordinates (analogous to CIELAB’s a/b axes) via `jab`  
//! - **Perceptually uniform “Ja′b′”** coordinates (CAM16-UCS) via `jab_prime`  
//! - **Inverse transform** back to XYZ with optional new white or viewing conditions via `xyz`.  
//! - **Chromatic adaptation** using CIECAT02 matrices  
//! - Utility functions for achromatic response, eccentricity, and more  

use super::{CamJCh, CamTransforms, ViewConditions};

use nalgebra::Vector3;

use crate::{
    error::Error,
    observer::Observer,
    rgb::{RgbSpace, WideRgb},
    xyz::{RelXYZ, XYZ},
};

/// CIECAM02 Color Appearance Model
///
/// Implements the CIECAM02 appearance model (CIE 248:2022), which predicts how we
/// perceive color under different viewing conditions. Building on CIECAM02 with
/// several enhancements, it computes the JCh correlates:
/// - **J**: Lightness  
/// - **C**: Chroma (colorfulness)  
/// - **h**: Hue angle  
///
/// It also performs chromatic adaptation, allowing colors to be accurately
/// transformed between different illuminants and viewing environments.
///
#[derive(Debug)]
pub struct CieCam02(CamJCh);

impl CieCam02 {
    /// Construct a CIECAM02 instance from precomputed JCh appearance correlates.
    ///
    /// This skips the usual XYZ→JCh conversion and initializes the model directly with:
    /// - J (lightness, 0…100)  
    /// - C (chroma, ≥ 0; in practice 0…~100+)  
    /// - h (hue angle in degrees, 0°…360°)  
    ///  
    /// # Arguments
    /// - `jch` — `[J, C, h]` correlates  
    /// - `xyzn` — Reference‐white `XYZ` (must share the same `Observer`)  
    /// - `vc` — Viewing conditions (`ViewConditions`)  
    ///
    /// # Returns
    /// A `CieCam02` instance initialized with the provided appearance correlates, white point, and viewing conditions.
    ///
    /// # Notes
    /// - Not every (J, C, h) triple corresponds to a real color.  
    /// - This constructor does _not_ validate the inputs; to ensure validity, you would need to  
    ///   perform the inverse transform (`.xyz()`) or convert to RGB and check for out-of-gamut values.  
    pub fn new(jch: [f64; 3], xyzn: XYZ, vc: ViewConditions) -> Self {
        Self(CamJCh {
            jch: Vector3::from(jch),
            xyzn,
            vc,
        })
    }

    /// Creates a new CieCam02 instance from the given XYZ tristimulus values of the color and the reference white,
    /// and the viewing conditions.
    /// The XYZ values must be share the same observer, otherwise an error is returned.
    /// # Arguments
    /// * `xyz` - The XYZ tristimulus values of the color to be transformed.
    /// * `xyzn` - The XYZ tristimulus values of the reference white.
    /// * `vc` - The viewing conditions to be used for the transformation.
    /// # Returns
    /// A Result containing the CieCam02 instance if successful, or a CmtError if an error occurs.
    /// # Errors
    /// Returns an error if the XYZ values are not in the same observer system.
    pub fn from_xyz(rxyz: RelXYZ, vc: ViewConditions) -> Self {
        Self(CamJCh::from_xyz(rxyz, vc, super::Cam::CieCam02))
    }

    /// Inverse-transform CIECAM02 appearance correlates back to CIE XYZ tristimulus values.
    ///
    /// You can optionally supply a new reference white or alternate viewing conditions:
    /// - `white_opt`:  
    ///   - `Some(white_xyz)`: use this `XYZ` as the white point for chromatic adaptation.  
    ///   - `None`: reuse the original white reference from when this `CieCam02` was created.
    /// - `vc_opt`:  
    ///   - `Some(new_vc)`: apply these `ViewConditions` instead of the original.  
    ///   - `None`: reuse the original viewing conditions.
    ///
    /// This is useful for:  
    /// - Matching colors across different illuminants or display environments.  
    /// - Verifying the round-trip accuracy of the forward and inverse transforms.
    ///
    /// # Arguments
    /// - `white_opt: Option<XYZ>`  
    /// - `vc_opt: Option<ViewConditions>`  
    ///
    /// # Returns
    /// - `Ok(XYZ)`: the reconstructed tristimulus values under the specified (or original) conditions.  
    /// - `Err(CmtError::RequireSameObserver)`: if `white_opt` has a different `Observer` than `self`.
    pub fn xyz(
        &self,
        opt_xyzn: Option<XYZ>,
        opt_viewconditions: Option<ViewConditions>,
    ) -> Result<RelXYZ, Error> {
        self.0
            .xyz(opt_xyzn, opt_viewconditions, super::Cam::CieCam02)
    }

    /// Convert this CieCam02 instance to RGB in the specified color space.
    pub fn rgb(
        &self,
        rgbspace: RgbSpace,
        opt_viewconditions: Option<ViewConditions>,
    ) -> Result<WideRgb, Error> {
        self.0
            .rgb(rgbspace, opt_viewconditions, super::Cam::CieCam02)
    }
}

impl CamTransforms for CieCam02 {
    /// Returns the JCh appearance correlates of this CieCam02 instance.
    fn jch_vec(&self) -> &Vector3<f64> {
        &self.0.jch
    }

    /// Returns the viewing conditions of this CieCam02 instance.
    fn view_conditions(&self) -> &ViewConditions {
        &self.0.vc
    }

    /// Returns the observer of this CieCam02 instance.
    fn observer(&self) -> Observer {
        self.0.xyzn.observer
    }

    fn xyzn(&self) -> &Vector3<f64> {
        &self.0.xyzn.xyz
    }
}

#[cfg(test)]
mod cam02_test {
    use approx::assert_abs_diff_eq;
    use nalgebra::Matrix3;

    use crate::cam::{CamTransforms, CieCam02, ViewConditions, M16, M16INV};
    use crate::observer::Observer;
    use crate::xyz::{RelXYZ, XYZ};

    #[test]
    fn test_m16() {
        approx::assert_abs_diff_eq!(M16INV * M16, Matrix3::identity(), epsilon = 1E-8);
    }

    #[test]
    // Worked example from CIECAM02 documentation, CIE159:2004, p. 11
    fn test_worked_example() {
        // forward transform (XYZ -> JCh)
        let white_point = XYZ::new([98.88, 90.0, 32.03], Observer::Cie1931);
        let rxyz = RelXYZ::new([19.31, 23.93, 10.14], white_point);
        // Table 4, column 2, CIE 159:2004.
        // La = 20 cd/m2;
        let vc = ViewConditions::new(20.0, 18.0, 0.69, 1.0, 1.0, None);
        let cam = CieCam02::from_xyz(rxyz, vc);
        let jch = cam.jch_vec();
        let &[j, c, h] = jch.as_ref();
        assert_abs_diff_eq!(j, 47.6856, epsilon = 1E-4); // J
        assert_abs_diff_eq!(c, 36.0527, epsilon = 1E-4); // C
        assert_abs_diff_eq!(h, 185.3445, epsilon = 1E-4); // h

        // Table 4, column 3, CIE 159:2004.
        // La = 200 cd/m2;
        let vc = ViewConditions::new(200.0, 18.0, 0.69, 1.0, 1.0, None);
        let cam = CieCam02::from_xyz(rxyz, vc);
        let jch = cam.jch_vec();
        let &[j, c, h] = jch.as_ref();
        assert_abs_diff_eq!(j, 48.0314, epsilon = 1E-4); // J
        assert_abs_diff_eq!(c, 38.7789, epsilon = 1E-4); // C
        assert_abs_diff_eq!(h, 191.0452, epsilon = 1E-4); // h
    }
}

#[cfg(test)]
mod cam02_round_trip_tests {
    use crate::{
        cam::{CamTransforms, CieCam02, ViewConditions},
        observer::Observer::{self, Cie1931},
        xyz::{RelXYZ, XYZ},
    };
    use approx::assert_abs_diff_eq;

    #[test]
    fn xyz_jch_xyz_round_trip() {
        // pick a handful of representative XYZs
        let samples = &[
            [19.01, 20.00, 21.78],
            [41.24, 21.26, 1.93],
            [95.05, 100.00, 108.88],
            [50.00, 50.00, 50.0],
            [0.00, 0.00, 0.00],
            [1.00, 1.00, 1.00],
            [70.00, 50.00, 30.00],
            [30.00, 60.00, 90.00],
            [12.14, 28.56, 5.00],
            [5.00, 12.14, 28.56],
            [100.00, 100.00, 100.00],
            [50.00, 50.00, 50.00],
            [20.00, 30.00, 40.00],
            [10.00, 20.00, 30.00],
        ];

        for &xyz_arr in samples {
            // forward transform (XYZ -> JCh)
            let xyz = XYZ::new(xyz_arr, Observer::Cie1931);
            let xyz_d65 = Cie1931.xyz_d65();
            let rxyz = RelXYZ::from_xyz(xyz, xyz_d65).unwrap();
            let cam = CieCam02::from_xyz(rxyz,ViewConditions::default());
            let jch = cam.jch();

            // inverse (JCh -> XYZ)
            let cam_back = CieCam02::new(jch, xyz_d65, ViewConditions::default());
            let xyz_back = cam_back.xyz(None, None).unwrap();

            // compare original vs. round-tripped XYZ
            assert_abs_diff_eq!(rxyz, xyz_back, epsilon = 1e-6);
        }
    }
}
