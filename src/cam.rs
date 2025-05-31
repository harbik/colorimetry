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
//! use colorimetry::cam::CamTransforms;
//!
//! let sample = XYZ::new([60.7, 49.6, 10.3], Observer::Std1931);
//! let white  = XYZ::new([96.46, 100.0, 108.62], Observer::Std1931);
//! let vc     = ViewConditions::new(16.0, 1.0, 1.0, 0.69, 40.0, None);
//!
//! let cam = CieCam16::from_xyz(sample, white, vc).unwrap();
//! let jch = cam.jch();
//! println!("JCh: {:?}", jch);
//! ```

mod viewconditions;
pub use viewconditions::ViewConditions;

mod cam16;
pub use crate::cam::cam16::CieCam16;

use crate::observer::Observer;
use nalgebra::{matrix, SMatrix, Vector3};
use std::f64::consts::PI;

#[derive(Debug)]
pub struct CamJCh {
    /// Colorimetric Observer used.
    /// The standard requires use of the CIE1931 standard observer.
    observer: Observer,

    /// Correlates of Lightness, Chroma, and hue-angle
    jch: Vector3<f64>,

    /// Tristimulus values of the reference white beng use
    xyzn: Vector3<f64>,

    /// Viewing Conditions
    vc: ViewConditions,
}

/// Values used in CieCam forward and backward transformations, only dependent on the view
/// conditions and the reference white.
/// For a set of colors, viewed under the same conditions, these have to be calculated only once.
#[derive(Clone, Copy, Debug)]
pub struct ReferenceValues {
    n: f64,
    z: f64,
    nbb: f64,
    ncb: f64,
    d_rgb: [f64; 3],
    aw: f64,
    qu: f64, // see lum_adapt
}

pub trait CamTransforms {
    const P1C: f64 = 50_000.0 / 13.0;
    const P3: f64 = 21.0 / 20.0;
    const NOM: f64 = 1.0; // is listed in the standard as (2.0+Self::P3)*460.0/C16_3; // this is 1!! But his is in Luo step 3 as part of nominators
    const DEN1: f64;
    const DEN2: f64;
    const RCPR_9: f64 = 1.0 / 0.9;
    const UCS_KL: f64 = 1.0;
    const UCS_C1: f64 = 0.007;
    const UCS_C2: f64 = 0.0228;

    fn jch_vec(&self) -> &Vector3<f64>;
    fn view_conditions(&self) -> &ViewConditions;
    fn observer(&self) -> Observer;
    fn xyzn(&self) -> &Vector3<f64>;
    fn xyz2cam_rgb(xyz_vec: Vector3<f64>) -> Vector3<f64>;

    fn reference_values(&self) -> ReferenceValues {
        let vc = self.view_conditions();
        let mut rgb_w = Self::xyz2cam_rgb(*self.xyzn());
        let vcd = vc.dd();
        let yw = self.xyzn()[1];
        let d_rgb = rgb_w.map(|v| vcd * yw / v + 1.0 - vcd);
        let n = vc.yb / yw;
        let z = n.sqrt() + 1.48;
        let nbb = 0.725 * n.powf(-0.2);
        let ncb = nbb;
        rgb_w.component_mul_assign(&Vector3::from(d_rgb));
        let qu = 150f64.max(rgb_w[0].max(rgb_w[1]).max(rgb_w[2]));

        // rgb_paw
        rgb_w.apply(|q| vc.lum_adapt(q, 0.26, qu));
        //        println!("***RGBaw {rgb_w}");
        let aw = achromatic_rsp(rgb_w, nbb);
        //        println!("***aw {aw} qu {qu}");

        ReferenceValues {
            n,
            z,
            nbb,
            ncb,
            d_rgb: d_rgb.into(),
            aw,
            qu,
        }
    }

    /// Returns the JCh values of the color as an array.
    ///
    /// The JCh values are a lightness, chroma, and hue angle representation of the color.
    fn jch(&self) -> [f64; 3] {
        (*self.jch_vec()).into()
    }

    /// CIECAM16 “raw Jab” coordinates, wich are analogous to CIELAB’s *a*/*b* values.
    ///
    /// **Note:**  
    /// - These raw Jab coordinates are _not_ perceptually uniform.  
    /// - Euclidean distances in this space do **not** reliably correspond to visual color differences.  
    ///
    /// For a perceptually uniform version, use `jab_prime()`, which applies the CAM16-UCS non-linear  
    /// stretching (with constants `C1 = 0.007` and `C2 = 0.0228`) to produce `(J′, a′, b′)`.
    fn jab(&self) -> [f64; 3] {
        let &[jj, cc, h] = self.jch_vec().as_ref();
        let vc = self.view_conditions();
        let m = cc * vc.f_l().powf(0.25);
        let mprime = 1.0 / Self::UCS_C2 * (1.0 + Self::UCS_C2 * m).ln();
        [
            (1.0 + 100.0 * Self::UCS_C1) * jj / (1.0 + Self::UCS_C1 * jj),
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
    fn jab_prime(&self) -> [f64; 3] {
        let &[jj, cc, h] = self.jch_vec().as_ref();
        let vc = self.view_conditions();
        let m = cc * vc.f_l().powf(0.25);
        let mprime = 1.0 / Self::UCS_C2 * (1.0 + Self::UCS_C2 * m).ln();
        [
            (1.0 + 100.0 * Self::UCS_C1) * jj / (1.0 + Self::UCS_C1 * jj),
            mprime * (h * PI / 180.0).cos(),
            mprime * (h * PI / 180.0).sin(),
        ]
    }

    /// Calculates the CIECAM-UCS ΔE′ (prime) color difference between two colors.
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
    fn de_cam(&self, other: &Self) -> Result<f64, crate::Error> {
        if self.observer() != other.observer() {
            return Err(crate::Error::RequireSameObserver);
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
        if (Self::UCS_KL - 1.0).abs() < f64::EPSILON {
            super::math::distance(jabp1, jabp2)
        } else {
            todo!()
        }
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

#[allow(dead_code)]
const MCAT02: SMatrix<f64, 3, 3> = matrix![
     0.7328,  0.4296,  -0.1624;
    -0.7036,  1.6975,   0.0061;
     0.0030,  0.0136,   0.9834;
];

/// Inverse CIECAT02 Chromatic Adaptation as a Matrix
#[allow(dead_code)]
const MCAT02INV: SMatrix<f64, 3, 3> = matrix![
    1.096123820835514, 		-0.2788690002182872, 	0.18274517938277304;
    0.45436904197535916,	 0.4735331543074117,	0.0720978037172291;
    -0.009627608738429353, 	-0.005698031216113419,	1.0153256399545427;
];

#[allow(dead_code)]
const MHPE: SMatrix<f64, 3, 3> = matrix![
     0.38971, 0.68898, -0.07868;
    -0.22981, 1.18340,  0.04641;
     0.00000, 0.00000,  1.00000;
];

#[allow(dead_code)]
const MHPEINVLUO: SMatrix<f64, 3, 3> = matrix![
    1.910197, -1.112124,  0.201908;
    0.370950,  0.629054, -0.000008;
    0.000000,  0.000000,  1.000000;
];

#[allow(dead_code)]
const MHPEINV: SMatrix<f64, 3, 3> = matrix![
    1.9101968340520348, -1.1121238927878747,  0.20190795676749937;
    0.3709500882486886,  0.6290542573926132, -0.000008055142184359149;
    0.0,  				 0.0,  				  1.0;
];

#[allow(dead_code)]
const MCAT02INVLUO: SMatrix<f64, 3, 3> = matrix![
     1.096124, -0.278869, 0.182745;
     0.454369,  0.473533, 0.072098;
    -0.009628, -0.005698, 1.015326;
];
