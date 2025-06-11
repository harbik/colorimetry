//! # CIECAM16 Color Appearance Model
//!
//! This module implements the **CIECAM16** appearance model (CIE 248:2022).  
//! It converts CIE XYZ tristimulus values into perceptual correlates of lightness (J),  
//! chroma (C), and hue angle (h) under specified viewing conditions, and supports:
//! - **Raw “Jab”** coordinates (analogous to CIELAB’s a/b axes) via `jab`  
//! - **Perceptually uniform “Ja′b′”** coordinates (CAM16-UCS) via `jab_prime`  
//! - **Inverse transform** back to XYZ with optional new white or viewing conditions via `xyz`.  
//! - **Chromatic adaptation** using CIECAT02 matrices  
//! - Utility functions for achromatic response, eccentricity, and more  
//!
use super::{CamJCh, CamTransforms, ViewConditions};

use nalgebra::Vector3;

use crate::{
    error::Error,
    observer::Observer,
    rgb::{RgbSpace, WideRgb},
    xyz::{RelXYZ, XYZ},
};

/// CIECAM16 Color Appearance Model
///
/// Implements the CIECAM16 appearance model (CIE 248:2022), which predicts how we
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
pub struct CieCam16(CamJCh);

impl CieCam16 {
    /// Construct a CIECAM16 instance from precomputed JCh appearance correlates.
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
    /// A `CieCam16` instance initialized with the provided appearance correlates, white point, and viewing conditions.
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

    /// Creates a new CieCam16 instance from the given XYZ tristimulus values of the color and the reference white,
    /// and the viewing conditions.
    /// The XYZ values must be share the same observer, otherwise an error is returned.
    /// # Arguments
    /// * `xyz` - The XYZ tristimulus values of the color to be transformed.
    /// * `xyzn` - The XYZ tristimulus values of the reference white.
    /// * `vc` - The viewing conditions to be used for the transformation.
    /// # Returns
    /// A Result containing the CieCam16 instance if successful, or a CmtError if an error occurs.
    /// # Errors
    /// Returns an error if the XYZ values are not in the same observer system.
    pub fn from_xyz(rxyz: RelXYZ, vc: ViewConditions) -> Self {
        Self(CamJCh::from_xyz(rxyz, vc, super::Cam::CieCam16))
    }

    /// Inverse-transform CIECAM16 appearance correlates back to CIE XYZ tristimulus values.
    ///
    /// You can optionally supply a new reference white or alternate viewing conditions:
    /// - `white_opt`:  
    ///   - `Some(white_xyz)`: use this `XYZ` as the white point for chromatic adaptation.  
    ///   - `None`: reuse the original white reference from when this `CieCam16` was created.
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
            .xyz(opt_xyzn, opt_viewconditions, super::Cam::CieCam16)
    }

    /// Converts the CieCam16 instance to RGB in the specified RGB color space.
    pub fn rgb(
        &self,
        rgbspace: RgbSpace,
        opt_viewconditions: Option<ViewConditions>,
    ) -> Result<WideRgb, Error> {
        self.0
            .rgb(rgbspace, opt_viewconditions, super::Cam::CieCam16)
    }

    /// Calculates the CIECAM16-UCS ΔE′ (prime) color difference between two colors.
    ///
    /// This method converts each color to its CAM16-UCS Ja′b′ coordinates and then computes the
    /// Euclidean distance, and applies a non-linear transformation to closely matches perceived
    /// color differences.
    ///
    /// # Parameters
    /// - `other`: The `CieCam16` instance to compare against.
    ///
    /// # Returns
    /// The ΔE′ value as an `f64`, representing the perceptual color difference.
    ///
    /// # Errors
    /// Returns an error if the observers of the two colors do not match.
    ///
    /// # References
    /// - [Comprehensive color solutions: CAM16, CAT16, and CAM16-UCS](https://www.researchgate.net/publication/318152296_Comprehensive_color_solutions_CAM16_CAT16_and_CAM16-UCS/)
    ///
    pub fn de_ucs(&self, other: &Self) -> Result<f64, crate::Error> {
        if self.observer() != other.observer() {
            return Err(crate::Error::RequireSameObserver);
        }

        let jabp1 = self.jab_prime();
        let jabp2 = other.jab_prime();
        let de_ucs_prime = Self::delta_e_prime_from_jabp(jabp1.as_ref(), jabp2.as_ref());
        Ok(1.41 * de_ucs_prime.powf(0.63))
    }
}

impl CamTransforms for CieCam16 {
    /// Returns the JCh appearance correlates of this CieCam16 instance.
    fn jch_vec(&self) -> &Vector3<f64> {
        &self.0.jch
    }

    /// Returns the viewing conditions of this CieCam16 instance.
    fn view_conditions(&self) -> &ViewConditions {
        &self.0.vc
    }

    /// Returns the observer of this CieCam16 instance.
    fn observer(&self) -> Observer {
        self.0.xyzn.observer
    }

    fn xyzn(&self) -> &Vector3<f64> {
        &self.0.xyzn.xyz
    }
}

#[cfg(test)]
mod cam16_test {
    use approx::assert_abs_diff_eq;
    use nalgebra::Matrix3;

    use crate::cam::{CamTransforms, CieCam16, ViewConditions, M16, M16INV};
    use crate::observer::Observer;
    use crate::xyz::{RelXYZ, XYZ};

    #[test]
    fn test_m16() {
        approx::assert_abs_diff_eq!(M16INV * M16, Matrix3::identity(), epsilon = 1E-8);
    }

    #[test]
    fn test_worked_example() {
        // see section 7 CIE 248:2022
        let xyz = XYZ::new([60.70, 49.60, 10.29], Observer::Cie1931);
        let xyzn = XYZ::new([96.46, 100.0, 108.62], Observer::Cie1931);
        let rxyz = RelXYZ::new(xyz, xyzn).unwrap();
        let vc = ViewConditions::new(40.0, 16.0, 0.69, 1.0, 1.0, None);
        let cam = CieCam16::from_xyz(rxyz, vc);
        let &[j, c, h] = cam.jch_vec().as_ref();
        //println!("J:\t{j:?}\nC:\t{c:?}\nh:\t{h:?}");
        approx::assert_abs_diff_eq!(j, 70.4406, epsilon = 1E-4);
        approx::assert_abs_diff_eq!(c, 58.6035, epsilon = 1E-4);
        approx::assert_abs_diff_eq!(h, 57.9145, epsilon = 1E-4);

        // inverse transformation, with no change in white adaptation of viewing conditions.
        let xyz_rev = cam.xyz(None, Some(vc)).unwrap();
        assert_abs_diff_eq!(rxyz, xyz_rev, epsilon = 1E-4);
    }
}

#[cfg(test)]
mod cam16_round_trip_tests {
    use crate::observer::Observer::Cie1931;
    use crate::prelude::*;
    use crate::xyz::RelXYZ;
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
            let xyz = XYZ::new(xyz_arr, Cie1931);
            let xyz_d65 = Cie1931.xyz_d65();
            let rxyz = RelXYZ::new(xyz, xyz_d65).unwrap();
            let cam = CieCam16::from_xyz(rxyz, ViewConditions::default());
            let jch = cam.jch();

            // inverse (JCh -> XYZ)
            let cam_back = CieCam16::new(jch, xyz_d65, ViewConditions::default());
            let rxyz_back = cam_back.xyz(None, None).unwrap();

            // compare original vs. round-tripped XYZ
            assert_abs_diff_eq!(rxyz, rxyz_back, epsilon = 1e-6);
        }
    }
}

#[cfg(test)]
mod rgb_test {
    #[test]
    #[cfg(all(
        feature = "supplemental-observers",
        feature = "munsell",
        feature = "cie-illuminants"
    ))]
    fn rgb_match() {
        use crate::{
            cam::{ViewConditions, CIE248_HOME_SCREEN},
            colorant::Munsell,
            illuminant::LED_B2,
            observer::Observer::{Cie1931, Cie2015_10},
            rgb::RgbSpace::SRGB,
        };

        let paint = Munsell::try_new("5BG5/8").unwrap();
        let vc = ViewConditions::average_surround(6.0);
        let cam_paint = Cie2015_10.ciecam16(&LED_B2, &paint, vc);
        let rgb_2015 = cam_paint
            .rgb(SRGB, Some(CIE248_HOME_SCREEN))
            .unwrap()
            .compress();

        // Use spectral representation of the Cie2015_10 RGB pixel
        // to convert to CIE1931
        let xyz_1931 = Cie1931.xyz(&rgb_2015, None);
        let rgb_1931 = xyz_1931.rgb(SRGB).compress();
        let [r, g, b]: [u8; 3] = rgb_1931.into();
        assert!(r == 0 && g == 113 && b == 138);
    }
}
