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

mod cam16;
pub use crate::cam::cam16::CieCam16;

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
use nalgebra::{matrix, Matrix3, SMatrix, Vector3};

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
    use approx::assert_abs_diff_eq;

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
