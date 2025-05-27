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
//! ## Example
//! ```rust
//! use colorimetry::cam::CieCam16;
//! use colorimetry::cam::ViewConditions;
//! use colorimetry::xyz::XYZ;
//! use colorimetry::observer::Observer;
//!
//! let sample = XYZ::new([60.7, 49.6, 10.3], Observer::Std1931);
//! let white  = XYZ::new([96.46, 100.0, 108.62], Observer::Std1931);
//! let vc     = ViewConditions::new(16.0, 1.0, 1.0, 0.69, 40.0, None);
//!
//! let cam = CieCam16::from_xyz(sample, white, vc).unwrap();
//! let jch = cam.jch();
//! println!("JCh: {:?}", jch);
//! ```
//!
//! *Methods and internals marked `pub(crate)` have been omitted for brevity.*

mod viewconditions;
pub use viewconditions::ViewConditions;

use std::f64::consts::PI;

const P1C: f64 = 50_000.0 / 13.0;
const P3: f64 = 21.0 / 20.0;
const C16_3: f64 = 1_403.0;
const NOM: f64 = 1.0; // is listed in the standard as (2.0+Self::P3)*460.0/C16_3; // this is 1!! But his is in Luo step 3 as part of nominators
const DEN1: f64 = ((2.0 + P3) * 220.0) / C16_3;
const DEN2: f64 = (P3 * 6300.0 - 27.0) / C16_3;
const RCPR_9: f64 = 1.0 / 0.9;
const UCS_KL: f64 = 1.0;
const UCS_C1: f64 = 0.007;
const UCS_C2: f64 = 0.0228;
use nalgebra::{matrix, vector, Matrix3, SMatrix, Vector3};

use crate::{
    cam::viewconditions::ReferenceValues,
    error::CmtError,
    geometry::distance,
    prelude::Observer,
    traits::{Filter, Light},
    xyz::XYZ,
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
pub struct CieCam16 {
    /// Colorimetric Observer used.
    /// The standard requires use of the CIE1931 standard observer.
    pub(crate) observer: Observer,
    /// Correlates of Lightness, Chroma, and hue-angle
    pub(crate) jch: Vector3<f64>,
    /// Tristimulus values of the reference white beng use
    pub(crate) xyzn: Vector3<f64>,
    /// Viewing Conditions
    pub(crate) vc: ViewConditions, // 6 f64
}

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
        Self {
            jch: Vector3::from(jch),
            xyzn: xyzn.xyz,
            vc,
            observer: xyzn.observer,
        }
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
    pub fn from_xyz(xyz: XYZ, xyzn: XYZ, vc: ViewConditions) -> Result<Self, CmtError> {
        let xyz_vec = xyz.xyz;
        let xyzn_vec = xyzn.xyz;
        if xyz.observer != xyzn.observer {
            return Err(CmtError::RequireSameObserver);
        }
        let ReferenceValues {
            n,
            z,
            nbb,
            ncb,
            d_rgb,
            aw,
            qu,
        } = ReferenceValues::new(xyzn_vec, vc);
        let vcdd = vc.dd();
        let vcfl = vc.f_l();
        let mut rgb = M16 * xyz_vec;
        rgb.component_mul_assign(&Vector3::from(d_rgb));
        rgb.apply(|v| vc.lum_adapt(v, 0.26, qu));

        let ca = rgb[0] - 12.0 / 11.0 * rgb[1] + rgb[2] / 11.0;
        let cb = (rgb[0] + rgb[1] - 2. * rgb[2]) / 9.0;

        // calculate h in radians
        let mut h = cb.atan2(ca);
        if h < 0.0 {
            h += 2.0 * PI;
        }

        // calculate J = jj
        let jj = 100.0 * (achromatic_rsp(rgb, nbb) / aw).powf(vc.c * z);

        // calculate C = cc
        let et = 0.25f64 * ((h + 2.0).cos() + 3.8);
        let t = (50000.0 / 13.0 * ncb * vc.nc * et * (ca * ca + cb * cb).sqrt())
            / (rgb[0] + rgb[1] + 21.0 / 20.0 * rgb[2]);
        let cc = t.powf(0.9) * (jj / 100.).sqrt() * (1.64 - (0.29f64).powf(n)).powf(0.73);

        Ok(Self {
            vc,
            jch: Vector3::new(jj, cc, h * 180.0 / PI),
            observer: xyz.observer,
            xyzn: xyzn_vec,
        })
    }

    /// Returns the JCh values of the color as a `Vector3<f64>`.
    ///
    /// The JCh values are a lightness, chroma, and hue angle representation of the color.
    pub fn jch(&self) -> [f64; 3] {
        self.jch.into()
    }

    /// CIECAM16 “raw Jab” coordinates, wich are analogous to CIELAB’s *a*/*b* values.
    ///
    /// **Note:**  
    /// - These raw Jab coordinates are _not_ perceptually uniform.  
    /// - Euclidean distances in this space do **not** reliably correspond to visual color differences.  
    ///
    /// For a perceptually uniform version, use `jab_prime()`, which applies the CAM16-UCS non-linear  
    /// stretching (with constants `C1 = 0.007` and `C2 = 0.0228`) to produce `(J′, a′, b′)`.
    pub fn jab(&self) -> [f64; 3] {
        let &[jj, cc, h] = self.jch.as_ref();
        let m = cc * self.vc.f_l().powf(0.25);
        let mprime = 1.0 / UCS_C2 * (1.0 + UCS_C2 * m).ln();
        [
            (1.0 + 100.0 * UCS_C1) * jj / (1.0 + UCS_C1 * jj),
            mprime * (h * PI / 180.0).cos(),
            mprime * (h * PI / 180.0).sin(),
        ]
    }

    /// This function returns the CIECAM16 UCS "Ja'b'" values, which are a perceptually uniform representation  
    /// of color.
    ///  
    /// **Note:**  
    /// - **J′**: Lightness  
    /// - **a′**, **b′**: Chromatic components, derived from the non-linear CAM16-UCS transformation applied  
    ///   to the raw Jab values.
    ///
    /// Euclidean distance between colors in this space more reliably matches human-perceived color differences.  
    /// The constants `C1 = 0.007` and `C2 = 0.0228` are used for the transformation to UCS type.  
    ///  
    /// For the raw, non-uniform Jab values, see the `jab()` function.
    pub fn jab_prime(&self) -> [f64; 3] {
        let &[jj, cc, h] = self.jch.as_ref();
        let m = cc * self.vc.f_l().powf(0.25);
        let mprime = 1.0 / UCS_C2 * (1.0 + UCS_C2 * m).ln();
        [
            (1.0 + 100.0 * UCS_C1) * jj / (1.0 + UCS_C1 * jj),
            mprime * (h * PI / 180.0).cos(),
            mprime * (h * PI / 180.0).sin(),
        ]
    }

    /// Returns the JC'h' values of the color as an array.
    ///
    /// The JC'h' values are a lightness, chroma, and hue angle representation of the color,
    /// where:
    /// - **J** is the lightness (0 to 100)
    /// - **C'** is the chroma (0 or greater)
    /// - **h'** is the hue angle in degrees (0° to 360°)
    ///
    /// This method is similar to `jab_prime()`, but it uses the JCh representation instead of the raw Jab.
    /// It applies the CAM16-UCS non-linear stretching to produce a perceptually uniform representation.
    ///
    /// # Returns
    /// An array containing the JCh values: `[J, C', h']`.
    pub fn jch_prime(&self) -> [f64; 3] {
        let &[jj, cc, h] = self.jch.as_ref();
        let m = cc * self.vc.f_l().powf(0.25);
        let mprime = 1.0 / UCS_C2 * (1.0 + UCS_C2 * m).ln();
        [
            (1.0 + 100.0 * UCS_C1) * jj / (1.0 + UCS_C1 * jj),
            mprime * (h * PI / 180.0).cos(),
            mprime * (h * PI / 180.0).sin(),
        ]
    }

    /// Calculates the CIECAM16-UCS ΔE′ (prime) color difference between two colors.
    ///
    /// This method converts each color to its CAM16-UCS Ja′b′ coordinates and then
    /// computes the Euclidean distance, which closely matches perceived color differences.
    ///
    /// # Formula
    /// ```text
    /// ΔE′ = sqrt((J1′ - J2′)² + (a1′ - a2′)² + (b1′ - b2′)²)
    /// ```
    ///
    /// # Parameters
    /// - `other`: The `CieCam16` instance to compare against.
    ///
    /// # Returns
    /// The ΔE′ value as an `f64`, representing the perceptual color difference.
    ///
    /// # Errors
    /// Returns an error if the observers of the two colors do not match.
    pub fn ciede2016(&self, other: &Self) -> Result<f64, CmtError> {
        if self.observer != other.observer {
            return Err(CmtError::RequireSameObserver);
        }

        let jabp1 = self.jab_prime();
        let jabp2 = other.jab_prime();
        Ok(Self::delta_e_prime_from_jabp(
            jabp1.as_ref(),
            jabp2.as_ref(),
        ))
    }

    /// Calculates the CIECAM16 distance between two Ja'b' values.
    /// This is the CIECAM16 distance formula, which is used to calculate the
    /// perceptual difference between two colors in the CIECAM16 color space.
    fn delta_e_prime_from_jabp(jabp1: &[f64], jabp2: &[f64]) -> f64 {
        if (UCS_KL - 1.0).abs() < f64::EPSILON {
            distance(jabp1, jabp2)
        } else {
            todo!()
        }
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
    ///
    /// # Example
    /// ```rust
    /// use colorimetry::prelude::*;
    /// // Original CAM16 instance:
    /// let sample_xyz = XYZ::new([60.7, 49.6, 10.3], Observer::Std1931);
    /// let white_xyz  = XYZ::new([96.46, 100.0, 108.62], Observer::Std1931);
    /// let vc     = ViewConditions::new(16.0, 1.0, 1.0, 0.69, 40.0, None);
    /// let cam = CieCam16::from_xyz(sample_xyz, white_xyz, vc).unwrap();
    ///
    /// // Inverse under same conditions:
    /// let back_to_xyz = cam.xyz(None, None).unwrap();
    ///
    /// // Inverse with a new white point:
    /// let new_white = XYZ::new([95.0, 100.0, 108.0], Observer::Std1931);
    /// let adapted_xyz = cam.xyz(Some(new_white), None).unwrap();
    /// ```
    pub fn xyz(
        &self,
        white_opt: Option<XYZ>,
        vc_opt: Option<ViewConditions>,
    ) -> Result<XYZ, CmtError> {
        let vc = vc_opt.unwrap_or(self.vc);
        let xyzn = if let Some(white) = white_opt {
            if white.observer == self.observer {
                white.xyz
            } else {
                return Err(CmtError::RequireSameObserver);
            }
        } else {
            self.xyzn
        };
        let ReferenceValues {
            n,
            z,
            nbb,
            ncb,
            d_rgb,
            aw,
            qu,
        } = ReferenceValues::new(xyzn, vc);
        let d_rgb_vec = Vector3::from(d_rgb);
        let &[lightness, chroma, hue_angle] = self.jch.as_ref();
        let t = (chroma / ((lightness / 100.0).sqrt() * (1.64 - 0.29f64.powf(n)).powf(0.73)))
            .powf(RCPR_9);
        let p1 = (P1C * vc.nc * ncb * eccentricity(hue_angle)) / t; // NaN if t=0, but OK, as check on t==0.0 if used
        let p2 = achromatic_response_from_lightness(aw, vc.c, z, lightness) / nbb + 0.305;
        let (a, b) = match hue_angle.to_radians().sin_cos() {
            (_, _) if t.is_nan() || t == 0.0 => (0.0, 0.0),
            (hs, hc) if hs.abs() >= hc.abs() => {
                let b = p2 * NOM / (p1 / hs + DEN1 * hc / hs + DEN2);
                (b * hc / hs, b)
            }
            (hs, hc) => {
                let a = p2 * NOM / (p1 / hc + DEN1 + DEN2 * hs / hc);
                (a, a * hs / hc)
            }
        };

        // rgb_a
        let m = matrix![ 460.0, 451.0, 288.0; 460.0, -891.0, -261.0; 460.0, -220.0, -6_300.0; ]
            / 1_403.0;
        let rgb_p = (m * vector![p2, a, b]).map(|x| inv_cone_adaptation(vc.f_l(), x)); // Step 4 & 5
        let rgb = rgb_p.component_div(&d_rgb_vec);

        let xyz = M16INV * rgb;
        Ok(XYZ::from_vecs(xyz, self.observer))
    }
}

#[inline]
pub fn achromatic_rsp(rgb: Vector3<f64>, nbb: f64) -> f64 {
    (2.0 * rgb[0] + rgb[1] + rgb[2] / 20.0 - 0.305) * nbb
}

pub fn eccentricity(hue_angle: f64) -> f64 {
    0.25 * ((hue_angle.to_radians() + 2.0).cos() + 3.8)
}

pub fn achromatic_response_from_lightness(aw: f64, c: f64, z: f64, lightness: f64) -> f64 {
    aw * (lightness / 100.0).powf(1.0f64 / (c * z))
}

fn inv_cone_adaptation(f_l: f64, x: f64) -> f64 {
    let x = x - 0.1;
    let t = 27.13 * x.abs() / (400.0 - x.abs());
    x.signum() * ((100.0 * t.powf(1.0 / 0.42)) / f_l)
}

const MCAT02: SMatrix<f64, 3, 3> = matrix![
     0.7328,  0.4296,  -0.1624;
    -0.7036,  1.6975,   0.0061;
     0.0030,  0.0136,   0.9834;
];

/**
   Inverse CIECAT02 Chromatic Adaptation as a Matrix
*/
const MCAT02INV: SMatrix<f64, 3, 3> = matrix![
    1.096123820835514, 		-0.2788690002182872, 	0.18274517938277304;
    0.45436904197535916,	 0.4735331543074117,	0.0720978037172291;
    -0.009627608738429353, 	-0.005698031216113419,	1.0153256399545427;
];

const MHPE: SMatrix<f64, 3, 3> = matrix![
     0.38971, 0.68898, -0.07868;
    -0.22981, 1.18340,  0.04641;
     0.00000, 0.00000,  1.00000;
];

const MHPEINVLUO: SMatrix<f64, 3, 3> = matrix![
    1.910197, -1.112124,  0.201908;
    0.370950,  0.629054, -0.000008;
    0.000000,  0.000000,  1.000000;
];

const MHPEINV: SMatrix<f64, 3, 3> = matrix![
    1.9101968340520348, -1.1121238927878747,  0.20190795676749937;
    0.3709500882486886,  0.6290542573926132, -0.000008055142184359149;
    0.0,  				 0.0,  				  1.0;
];

const MCAT02INVLUO: SMatrix<f64, 3, 3> = matrix![
     1.096124, -0.278869, 0.182745;
     0.454369,  0.473533, 0.072098;
    -0.009628, -0.005698, 1.015326;
];

const M16: Matrix3<f64> = matrix![
    0.401288, 0.650173, -0.051461;
    -0.250268, 1.204414, 0.045854;
    -0.002079, 0.048952, 0.953127
];

const M16INV: Matrix3<f64> = matrix![
    1.86206786, -1.01125463, 0.14918677;
    0.38752654, 0.62144744, -0.00897398;
    -0.01584150, -0.03412294, 1.04996444
];

const MRGBAINV: Matrix3<f64> = matrix![
460.0/C16_3, 451.0/C16_3, 288.0/C16_3;
460.0/C16_3, -891.0/C16_3, -261.0/C16_3;
460.0/C16_3, -220.0/C16_3, -6_300.0/C16_3;
 ];

#[cfg(test)]
mod cam_test {
    use super::*;
    use crate::prelude::*;
    use approx::{assert_abs_diff_eq, assert_ulps_eq};

    #[test]
    fn test_m16() {
        approx::assert_abs_diff_eq!(M16INV * M16, Matrix3::identity(), epsilon = 1E-8);
    }

    #[test]
    fn test_worked_example() {
        // see section 7 CIE 248:2022
        let xyz = XYZ::new([60.70, 49.60, 10.29], Observer::Std1931);
        let xyzn = XYZ::new([96.46, 100.0, 108.62], Observer::Std1931);
        let vc = ViewConditions::new(16.0, 1.0, 1.0, 0.69, 40.0, None);
        let cam = CieCam16::from_xyz(xyz, xyzn, vc).unwrap();
        let &[j, c, h] = cam.jch.as_ref();
        // println!("J:\t{j:?}\nC:\t{c:?}\nh:\t{h:?}");
        approx::assert_abs_diff_eq!(j, 70.4406, epsilon = 1E-4);
        approx::assert_abs_diff_eq!(c, 58.6035, epsilon = 1E-4);
        approx::assert_abs_diff_eq!(h, 57.9145, epsilon = 1E-4);

        // inverse transformation, with no change in white adaptation of viewing conditions.
        let xyz_rev = cam.xyz(None, None).unwrap();
        assert_abs_diff_eq!(xyz, xyz_rev, epsilon = 1E-4);
    }
}

#[cfg(test)]
mod round_trip_tests {
    use crate::prelude::*;
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
            let xyz = XYZ::new(xyz_arr, Observer::Std1931);
            let xyz_d65 = CIE1931.xyz_d65();
            let cam = CieCam16::from_xyz(xyz, xyz_d65, ViewConditions::default()).unwrap();
            let jch = cam.jch();

            // inverse (JCh -> XYZ)
            let cam_back = CieCam16::new(jch, xyz_d65, ViewConditions::default());
            let xyz_back = cam_back.xyz(None, None).unwrap();

            // compare original vs. round-tripped XYZ
            let orig = XYZ::new(xyz_arr, Observer::Std1931);
            let [x0, y0, z0] = orig.values();
            let [x1, y1, z1] = xyz_back.values();

            assert_abs_diff_eq!(x0, x1, epsilon = 1e-6);
            assert_abs_diff_eq!(y0, y1, epsilon = 1e-6);
            assert_abs_diff_eq!(z0, z1, epsilon = 1e-6);
        }
    }
}
