//! # CIE L\*a\*b\* Color Model
//!
//! This module provides a **CIE L\*a\*b\*** representation (`CieLab`) and related
//! utilities for color difference metrics in Rust:
//!
//! - **`CieLab` struct**  
//!   Encapsulates L\*, a\*, b\* coordinates together with a reference white (XYZ) and
//!   a standard observer.  
//!
//! - **Construction**  
//!   - `CieLab::new([L, a, b], white_xyz)` — create directly from L\*a\*b\* values.  
//!   - `CieLab::from_xyz(xyz, white_xyz)` — convert from CIE XYZ under a given white point.  
//!
//! - **Accessors**  
//!   - `.values()` — returns `[L, a, b]`.  
//!   - AsRef<[f64; 3]> for ergonomic tuple-like access.  
//!
//! - **Color‐difference methods**  
//!   - `ciede(&other)` — plain Euclidean ΔE\*ab (common but not perceptually uniform).  
//!   - `ciede2000(&other)` — advanced CIEDE2000 formula (better matches human vision).  
//!
//! - **Error handling**  
//!   Ensures both colors share the same observer and illuminant, returning
//!   `CmtError::RequireSameObserver` or `CmtError::RequiresSameIlluminant` otherwise.
//!
//! - **WASM bindings**  
//!   Exported via `#[wasm_bindgen]` for JavaScript interoperability.
//!

use approx::{abs_diff_eq, ulps_eq, AbsDiffEq};
use nalgebra::Vector3;
use std::f64::consts::PI;

mod lch;
pub use lch::CieLCh;

mod gamut;
pub use gamut::CieLChGamut;

use crate::{
    error::Error,
    xyz::{RelXYZ, XYZ},
};

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CieLab {
    lab: Vector3<f64>,
    white_point: XYZ, // Reference white tristimulus value
}

impl CieLab {
    /// Creates a new CIE L*a*b* color from the given L*a*b* values and reference white.
    /// # Arguments
    /// * `lab` - The L*a*b* color values as an array of three f64 values.
    /// * `white_point` - The reference white tristimulus value.
    ///
    /// # Returns
    /// A new `CieLab` instance.
    pub fn new(lab: [f64; 3], xyzn: XYZ) -> CieLab {
        let lab = Vector3::from(lab);
        CieLab {
            lab,
            white_point: xyzn,
        }
    }

    pub fn a(&self) -> f64 {
        self.lab[1]
    }

    pub fn b(&self) -> f64 {
        self.lab[2]
    }

    pub fn l(&self) -> f64 {
        self.lab[0]
    }

    /// Creates a new CIE L*a*b* color from the given XYZ color and reference white.
    ///
    /// # Arguments
    /// * `rxyz` - The RelXYZ color to convert.
    ///
    /// # Returns
    /// A `Result` containing the CIE L*a*b* color or an error if the observers do not match.
    pub fn from_rxyz(xyz: RelXYZ) -> CieLab {
        CieLab {
            lab: lab(xyz.xyz().xyz, xyz.white_point().xyz),
            white_point: xyz.white_point(),
        }
    }

    /// Converts the CIE L\*a\*b\* color back to XYZ using the reference white.
    /// # Returns
    /// The XYZ color corresponding to the CIE L\*a\*b\* values.
    pub fn rxyz(&self) -> RelXYZ {
        let xyz = xyz_from_cielab(self.lab, self.white_point.xyz);
        RelXYZ::from_vec(xyz, self.white_point)
    }

    /// Converts the CIE L\*a\*b\* color back to XYZ using the reference white.
    /// # Returns
    /// The XYZ color corresponding to the CIE L\*a\*b\* values.
    pub fn xyz(&self) -> XYZ {
        self.rxyz().xyz()
    }

    /// Returns the reference white tristimulus value for this CIE L*a*b* color.
    pub fn white_point(&self) -> XYZ {
        self.white_point
    }

    /// Sets the reference white luminance for this CIE L*a*b* color, in units of cd/m².
    ///
    /// # Arguments
    /// * `luminance` - The desired luminance level for the reference white.
    /// # Returns
    /// A new `CieLab` instance with the adjusted luminance.     
    ///
    /// This adjusts the reference white to the specified illuminance level.
    ///
    /// # Notes
    /// - This does not change the L\*a\*b\* values directly; it modifies the reference white
    ///   to scale the white reference luminance.
    /// - Typically the value is set to 100 for normalized white luminance, but more advanced models,
    ///   such as CIECAM16, may use different luminance levels for perceptual accuracy.
    pub fn set_white_luminance(mut self, luminance: f64) -> CieLab {
        // Adjust the reference white to the new luminance
        let scale = luminance / self.white_point.xyz.y;
        self.white_point.xyz *= scale;
        self
    }

    /// Computes the Euclidean ΔE*ab color difference between two CIE L\*a\*b\* colors.
    ///
    /// This function measures the straight-line distance in L\*a\*b\* space:
    /// ΔE = sqrt((L₁−L₂)² + (a₁−a₂)² + (b₁−b₂)²)
    ///
    /// # Arguments
    /// - `other` – The second Lab color to compare against.
    ///
    /// # Returns
    /// - `Ok(de)` – The ΔE value if both colors share the same observer and illuminant.
    ///
    /// # Errors
    /// - `CmtError::RequireSameObserver` if the two colors use different standard observers.  
    /// - `CmtError::IlluminantMismatch` if they were computed under different illuminants.  
    ///
    /// # Notes
    /// The plain Euclidean ΔE*ab is commonly used but does not always match perceived differences
    /// as well as more advanced formulas (e.g., CIEDE2000, or CIECAM16DE).
    /// # Example
    /// ```
    /// use colorimetry::{lab::CieLab, xyz::{RelXYZ, XYZ}, observer::Observer, Error};
    ///
    /// let xyz1 = XYZ::new([36.0, 70.0, 12.0], Observer::Cie1931);
    /// let xyz2 = XYZ::new([35.0, 71.0, 11.0], Observer::Cie1931);
    /// let lab1 = CieLab::from_rxyz(RelXYZ::with_d65(xyz1));
    /// let lab2 = CieLab::from_rxyz(RelXYZ::with_d65(xyz2));
    /// let de = lab1.ciede(&lab2).unwrap();
    /// # approx::assert_abs_diff_eq!(de, 6.57, epsilon = 0.01);
    /// //  ΔE=6.57
    /// ```
    pub fn ciede(&self, other: &Self) -> Result<f64, Error> {
        if self.white_point.observer != other.white_point.observer {
            return Err(Error::RequireSameObserver);
        }
        if ulps_eq!(self.white_point, other.white_point) {
            let &[l1, a1, b1] = self.lab.as_ref();
            let &[l2, a2, b2] = other.lab.as_ref();
            Ok(((l2 - l1).powi(2) + (a2 - a1).powi(2) + (b2 - b1).powi(2)).sqrt())
        } else {
            Err(Error::RequiresSameIlluminant)
        }
    }

    ///
    /// Computes the CIEDE2000 ΔE color difference between two CIE L*a*b* colors.
    ///
    /// This is a more advanced formula that accounts for perceptual non-uniformities in the L*a*b* space.
    /// # Arguments
    /// - `other` – The second Lab color to compare against.
    ///
    /// # Returns
    /// - `Ok(de)` – The ΔE value if both colors share the same observer and illuminant.
    ///
    /// # Errors
    /// - `CmtError::RequireSameObserver` if the two colors use different standard observers.
    /// - `CmtError::RequiresSameIlluminant` if they were computed under different illuminants.
    ///
    /// # Notes
    /// CIEDE2000 is generally preferred over the plain Euclidean ΔE*ab for color difference calculations,
    /// as it better matches human perception of color differences.
    /// # Example
    /// ```
    /// use colorimetry::prelude::*;
    ///
    /// // Sharma et al. (2005) test case 25
    /// let xyz_d65 = Cie1931.xyz_d65();
    /// let lab1 = CieLab::new([60.2574, -34.0099, 36.2677], xyz_d65);
    /// let lab2 = CieLab::new([60.4626, -34.1751, 39.4387], xyz_d65);
    /// let de = lab1.ciede2000(&lab2).unwrap();
    /// approx::assert_abs_diff_eq!(de, 1.2644, epsilon = 1E-4);
    /// ```
    pub fn ciede2000(&self, other: &Self) -> Result<f64, Error> {
        if self.white_point.observer != other.white_point.observer {
            return Err(Error::RequireSameObserver);
        }
        if ulps_eq!(self.white_point, other.white_point) {
            Ok(delta_e_ciede2000(self.lab, other.lab))
        } else {
            Err(Error::RequiresSameIlluminant)
        }
    }

    /// Returns the CIE L*a*b* color values as an array of three f64 values.
    /// # Returns
    /// An array containing the L*, a*, and b* values of the color.
    pub fn values(&self) -> [f64; 3] {
        *self.lab.as_ref()
    }

    /// Validates the CIE L*a*b* color values.
    /// # Returns
    /// `true` if the L*, a*, and b* values are within valid ranges and
    /// the round-trip conversion to and from XYZ is consistent; `false` otherwise.
    pub fn is_valid(&self) -> bool {
        let [l, a, b] = self.values();
        if (0.0..=100.0).contains(&l)
            && (-200.0..=200.0).contains(&a)
            && ((-200.0..=200.0).contains(&b))
        {
            // check if the associated XYZ is valid
            self.rxyz().is_valid()
        } else {
            false
        }
    }

    pub fn is_black(&self, epsilon: f64) -> bool {
        let mut black = *self;
        black.lab = Vector3::new(0.0, 0.0, 0.0);
        let de = self.ciede2000(&black).unwrap(); // same observer and white point
        de < epsilon
    }
}

impl AsRef<[f64; 3]> for CieLab {
    fn as_ref(&self) -> &[f64; 3] {
        self.lab.as_ref()
    }
}

impl AbsDiffEq for CieLab {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-6
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.white_point.observer == other.white_point.observer
            && self.white_point == other.white_point
            && abs_diff_eq!(self.lab, other.lab, epsilon = epsilon)
    }
}

/// Convert CIEXYZ to CIELAB (D50, D65 or any other reference white passed in `xyzn`)
///
/// * `xyz`  – XYZ vector
/// * `xyzn` – reference-white Xn, Yn, Zn (must be in the **same scale** as the XYZ values you want
///   out – e.g. Yn = 100 for 0-to-100 “percent” data or Yn = 1.0 for ICC/PCS-scaled values)
/// # Returns
/// A vector containing the CIELAB L*, a*, b* values.
/// # References
/// - [CIE 15:2018](https://www.cie.co.at/publications/cie-15-2018-colorimetry), equation 8.3 to 8.11
fn lab(xyz: Vector3<f64>, xyzn: Vector3<f64>) -> Vector3<f64> {
    const DELTA: f64 = 0.008856451679035631; // 24f64 / 116.0).powi(3);
    const LABC1: f64 = 841.0 / 108.0;
    const LABC2: f64 = 16.0 / 116.0;

    let &[x, y, z] = xyz.as_ref();
    let &[xn, yn, zn] = xyzn.as_ref();

    let lab_f = |t: f64| {
        if t > DELTA {
            t.cbrt()
        } else {
            LABC1 * t + LABC2
        }
    };

    Vector3::new(
        116f64 * lab_f(y / yn) - 16f64,
        500f64 * (lab_f(x / xn) - lab_f(y / yn)),
        200f64 * (lab_f(y / yn) - lab_f(z / zn)),
    )
}

/// Convert CIE L*a*b* to CIEXYZ (D50, D65 or any other reference white passed in `xyzn`)
/// * `lab`  – CIE L*a*b* vector
/// * `xyzn` – reference-white Xn, Yn, Zn (must be in the **same scale** as the XYZ values you want
///   out – e.g. Yn = 100 for 0-to-100 “percent” data or Yn = 1.0 for ICC/PCS-scaled values)
/// # Returns
/// A vector containing the CIEXYZ X, Y, Z values.
/// /// # References
/// - [CIE 15:2018](https://www.cie.co.at/publications/colorimetry-4th-edition), equation 8.3 to 8.11
fn xyz_from_cielab(lab: Vector3<f64>, xyzn: Vector3<f64>) -> Vector3<f64> {
    // CIE constants
    const DELTA_INV: f64 = 24f64 / 116.0;
    const LABC1INV: f64 = 108.0 / 841.0; // LABC1
    const LABC2: f64 = 16.0 / 116.0;

    let &[l, a, b] = lab.as_ref();
    let &[xn, yn, zn] = xyzn.as_ref();

    let fy = (l + 16.0) / 116.0; // C.1
    let fx = a / 500.0 + fy; // C.2
    let fz = fy - b / 200.0; // C.3

    // C.4 - C.9
    let f_inv = |f: f64| {
        if f > DELTA_INV {
            // If the value is above the threshold (DELTA), use the cube of the value.
            // Otherwise, apply the linear transformation to approximate the non-linear behavior.
            f.powi(3)
        } else {
            (f - LABC2) * LABC1INV
        }
    };

    // C.4 - C.9
    let x = xn * f_inv(fx);
    let y = yn * f_inv(fy);
    let z = zn * f_inv(fz);

    Vector3::new(x, y, z)
}

/// Compute the CIEDE2000 ΔE between two CIE L*a*b* triples.
fn delta_e_ciede2000(lab1: Vector3<f64>, lab2: Vector3<f64>) -> f64 {
    // unpack
    let (l1, a1, b1) = (lab1[0], lab1[1], lab1[2]);
    let (l2, a2, b2) = (lab2[0], lab2[1], lab2[2]);

    // Step 1: C* and h*
    let c1 = (a1 * a1 + b1 * b1).sqrt();
    let c2 = (a2 * a2 + b2 * b2).sqrt();

    // Step 2: G factor
    let c_bar = (c1 + c2) / 2.0;
    let c_bar7 = c_bar.powi(7);
    let g = 0.5 * (1.0 - (c_bar7 / (c_bar7 + 25_f64.powi(7))).sqrt());

    // Step 3: a' and C'
    let a1p = a1 * (1.0 + g);
    let a2p = a2 * (1.0 + g);
    let c1p = (a1p * a1p + b1 * b1).sqrt();
    let c2p = (a2p * a2p + b2 * b2).sqrt();

    // Step 4: h'
    let h1p = b1.atan2(a1p).rem_euclid(2.0 * PI);
    let h2p = b2.atan2(a2p).rem_euclid(2.0 * PI);

    // Step 5: ΔL' and ΔC'
    let delta_lp = l2 - l1;
    let delta_cp = c2p - c1p;

    // Step 6: ΔH'
    let delta_hp = {
        let mut dh = h2p - h1p;
        if dh.abs() > PI {
            if dh > 0.0 {
                dh -= 2.0 * PI
            } else {
                dh += 2.0 * PI
            }
        }
        2.0 * (c1p * c2p).sqrt() * (dh / 2.0).sin()
    };

    // Step 7: Means
    let l_bar_p = (l1 + l2) / 2.0;
    let c_bar_p = (c1p + c2p) / 2.0;
    let h_bar_p = if (h1p - h2p).abs() > PI {
        // wrap‐around average
        let sum = h1p + h2p + 2.0 * PI;
        (sum / 2.0).rem_euclid(2.0 * PI)
    } else {
        (h1p + h2p) / 2.0
    };

    // Step 8: T parameter
    let t = 1.0 - 0.17 * (h_bar_p - 30.0_f64.to_radians()).cos()
        + 0.24 * (2.0 * h_bar_p).cos()
        + 0.32 * (3.0 * h_bar_p + 6.0_f64.to_radians()).cos()
        - 0.20 * (4.0 * h_bar_p - 63.0_f64.to_radians()).cos();

    // Step 9: Δθ, R_C, S_L/C/H, R_T
    let delta_theta =
        30.0_f64.to_radians() * (-((h_bar_p.to_degrees() - 275.0) / 25.0).powi(2)).exp();
    let rc = 2.0 * (c_bar_p.powi(7) / (c_bar_p.powi(7) + 25_f64.powi(7))).sqrt();
    let sl = 1.0 + (0.015 * (l_bar_p - 50.0).powi(2)) / (20.0 + (l_bar_p - 50.0).powi(2)).sqrt();
    let sc = 1.0 + 0.045 * c_bar_p;
    let sh = 1.0 + 0.015 * c_bar_p * t;
    let rt = -rc * (2.0 * delta_theta).sin();

    // Step 10: Final ΔE
    ((delta_lp / sl).powi(2)
        + (delta_cp / sc).powi(2)
        + (delta_hp / sh).powi(2)
        + rt * (delta_cp / sc) * (delta_hp / sh))
        .sqrt()
}

#[cfg(test)]
mod tests {
    use crate::{
        illuminant::CieIlluminant,
        lab::CieLab,
        observer::Observer::Cie1931,
        xyz::{RelXYZ, XYZ},
    };
    use approx::assert_abs_diff_eq;
    use nalgebra::vector;

    #[test]
    fn lab_roundtrip_test() {
        let xyz_values = [36.0, 70.0, 12.0];
        let xyz = XYZ::new(xyz_values, Cie1931);
        let rxyz = RelXYZ::with_d65(xyz);
        let lab = CieLab::from_rxyz(rxyz);
        let rxyz_back = lab.rxyz();
        assert_abs_diff_eq!(rxyz, rxyz_back, epsilon = 1e-10);
    }

    #[test]
    #[cfg(feature = "munsell")]
    fn test_lab_roundtrip_munsell() {
        use crate::colorant::MunsellCollection;
        MunsellCollection.into_iter().for_each(|colorant| {
            let xyz = Cie1931.rel_xyz(&CieIlluminant::D65, &colorant);
            let lab = CieLab::from_rxyz(xyz);
            let xyz_back = lab.rxyz();
            assert_abs_diff_eq!(xyz, xyz_back, epsilon = 1e-10);
        });
    }

    #[test]
    fn test_lab_roundtrip_monochromes() {
        let monos = Cie1931.monochromes(CieIlluminant::D65);
        monos.into_iter().for_each(|(_l, xyz)| {
            let lab = CieLab::from_rxyz(xyz);
            let xyz_back = lab.rxyz();
            assert_abs_diff_eq!(xyz, xyz_back, epsilon = 1e-10);
        });
    }

    #[test]
    fn delta_e_ciede2000_example1() {
        // Example 1 from Sharma et al. (2005):
        // ΔE₀₀ between (50.0000,  2.6772, –79.7751) and (50.0000,  0.0000, –82.7485) ≈ 2.0425
        let lab1 = vector![50.0000, 2.6772, -79.7751];
        let lab2 = vector![50.0000, 0.0000, -82.7485];
        let de = super::delta_e_ciede2000(lab1, lab2);
        assert_abs_diff_eq!(de, 2.0425, epsilon = 1e-4);
    }

    #[test]
    fn delta_e_ciede2000_example2() {
        // Example 2 from Sharma et al. (2005):
        // ΔE₀₀ between (50.0000,  3.1571, –77.2803) and (50.0000,  0.0000, –82.7485) ≈ 2.8615
        let lab1 = vector![50.0000, 3.1571, -77.2803];
        let lab2 = vector![50.0000, 0.0000, -82.7485];
        let de = super::delta_e_ciede2000(lab1, lab2);
        assert_abs_diff_eq!(de, 2.8615, epsilon = 1e-4);
    }

    #[test]
    fn supplemental_dataset() {
        // Sharma et al. (2005) supplemental test data: [L1, a1, b1], [L2, a2, b2], ΔE₀₀
        let cases = [
            ([50.0, 2.6772, -79.7751], [50.0, 0.0, -82.7485], 2.0425),
            ([50.0, 3.1571, -77.2803], [50.0, 0.0, -82.7485], 2.8615),
            ([50.0, 2.8361, -74.0200], [50.0, 0.0, -82.7485], 3.4412),
            ([50.0, -1.3802, -84.2814], [50.0, 0.0, -82.7485], 1.0000),
            ([50.0, -1.1848, -84.8006], [50.0, 0.0, -82.7485], 1.0000),
            ([50.0, -0.9009, -85.5211], [50.0, 0.0, -82.7485], 1.0000),
            ([50.0, 0.0000, 0.0000], [50.0, -1.0000, 2.0000], 2.3669),
            ([50.0, -1.0000, 2.0000], [50.0, 0.0000, 0.0000], 2.3669),
            ([50.0, 2.4900, -0.0010], [50.0, -2.4900, 0.0009], 7.1792),
            ([50.0, 2.4900, -0.0010], [50.0, -2.4900, 0.0010], 7.1792),
            ([50.0, 2.4900, -0.0010], [50.0, -2.4900, 0.0011], 7.2195),
            ([50.0, 2.4900, -0.0010], [50.0, -2.4900, 0.0012], 7.2195),
            ([50.0, -0.0010, 2.4900], [50.0, 0.0009, -2.4900], 4.8045),
            ([50.0, -0.0010, 2.4900], [50.0, 0.0010, -2.4900], 4.8045),
            ([50.0, -0.0010, 2.4900], [50.0, 0.0011, -2.4900], 4.7461),
            ([50.0, 2.5000, 0.0000], [50.0, 0.0000, -2.5000], 4.3065),
            ([50.0, 2.5000, 0.0000], [73.0, 25.0000, -18.0000], 27.1492),
            ([50.0, 2.5000, 0.0000], [61.0, -5.0000, 29.0000], 22.8977),
            ([50.0, 2.5000, 0.0000], [56.0, -27.0000, -3.0000], 31.9030),
            ([50.0, 2.5000, 0.0000], [58.0, 24.0000, 15.0000], 19.4535),
            ([50.0, 2.5000, 0.0000], [50.0, 3.1736, 0.5854], 1.0000),
            ([50.0, 2.5000, 0.0000], [50.0, 3.2972, 0.0000], 1.0000),
            ([50.0, 2.5000, 0.0000], [50.0, 1.8634, 0.5757], 1.0000),
            ([50.0, 2.5000, 0.0000], [50.0, 3.2592, 0.3350], 1.0000),
            (
                [60.2574, -34.0099, 36.2677],
                [60.4626, -34.1751, 39.4387],
                1.2644,
            ),
            (
                [63.0109, -31.0961, -5.8663],
                [62.8187, -29.7946, -4.0864],
                1.2630,
            ),
            (
                [61.2901, 3.7196, -5.3901],
                [61.4292, 2.2480, -4.9620],
                1.8731,
            ),
            (
                [35.0831, -44.1164, 3.7933],
                [35.0232, -40.0716, 1.5901],
                1.8645,
            ),
            (
                [22.7233, 20.0904, -46.6940],
                [23.0331, 14.9730, -42.5619],
                2.0373,
            ),
            (
                [36.4612, 47.8580, 18.3852],
                [36.2715, 50.5065, 21.2231],
                1.4146,
            ),
            (
                [90.8027, -2.0831, 1.4410],
                [91.1528, -1.6435, 0.0447],
                1.4441,
            ),
            (
                [90.9257, -0.5406, -0.9208],
                [88.6381, -0.8985, -0.7239],
                1.5381,
            ),
            (
                [6.7747, -0.2908, -2.4247],
                [5.8714, -0.0985, -2.2286],
                0.6377,
            ),
            (
                [2.0776, 0.0795, -1.1350],
                [0.9033, -0.0636, -0.5514],
                0.9082,
            ),
        ];

        for &(lab1_arr, lab2_arr, expected) in &cases {
            let lab1 = vector![lab1_arr[0], lab1_arr[1], lab1_arr[2]];
            let lab2 = vector![lab2_arr[0], lab2_arr[1], lab2_arr[2]];
            let de = super::delta_e_ciede2000(lab1, lab2);
            assert_abs_diff_eq!(de, expected, epsilon = 1e-4);
        }
    }

    #[test]
    fn supplemental_euclidean_delta_e() {
        // Supplemental dataset from Sharma et al. (2005), Table I:
        // ([L1, a1, b1], [L2, a2, b2], expected ΔE*ab)
        let cases: &[([f64; 3], [f64; 3], f64)] = &[
            ([50.0, 2.6772, -79.7751], [50.0, 0.0, -82.7485], 4.0011),
            ([50.0, 3.1571, -77.2803], [50.0, 0.0, -82.7485], 6.3142),
            ([50.0, 2.8361, -74.0200], [50.0, 0.0, -82.7485], 9.1777),
            ([50.0, -1.3802, -84.2814], [50.0, 0.0, -82.7485], 2.0627),
            ([50.0, -1.1848, -84.8006], [50.0, 0.0, -82.7485], 2.3696),
            ([50.0, -0.9009, -85.5211], [50.0, 0.0, -82.7485], 2.9153),
            ([50.0, 0.0000, 0.0000], [50.0, -1.0000, 2.0000], 2.2361),
            ([50.0, -1.0000, 2.0000], [50.0, 0.0000, 0.0000], 2.2361),
            ([50.0, 2.4900, -0.0010], [50.0, -2.4900, 0.0009], 4.9800),
            ([50.0, 2.4900, -0.0010], [50.0, -2.4900, 0.0010], 4.9800),
            ([50.0, 2.4900, -0.0010], [50.0, -2.4900, 0.0011], 4.9800),
            ([50.0, 2.4900, -0.0010], [50.0, -2.4900, 0.0012], 4.9800),
            ([50.0, -0.0010, 2.4900], [50.0, 0.0009, -2.4900], 4.9800),
            ([50.0, -0.0010, 2.4900], [50.0, 0.0010, -2.4900], 4.9800),
            ([50.0, -0.0010, 2.4900], [50.0, 0.0011, -2.4900], 4.9800),
            ([50.0, 2.5000, 0.0000], [50.0, 0.0000, -2.5000], 3.5355),
            ([50.0, 2.5000, 0.0000], [73.0, 25.0000, -18.0000], 36.8680),
            ([50.0, 2.5000, 0.0000], [61.0, -5.0000, 29.0000], 31.9100),
            ([50.0, 2.5000, 0.0000], [56.0, -27.0000, -3.0000], 30.2531),
            ([50.0, 2.5000, 0.0000], [58.0, 24.0000, 15.0000], 27.4089),
            ([50.0, 2.5000, 0.0000], [50.0, 3.1736, 0.5854], 0.8924),
            ([50.0, 2.5000, 0.0000], [50.0, 3.2972, 0.0000], 0.7972),
            ([50.0, 2.5000, 0.0000], [50.0, 1.8634, 0.5757], 0.8583),
            ([50.0, 2.5000, 0.0000], [50.0, 3.2592, 0.3350], 0.8298),
            (
                [60.2574, -34.0099, 36.2677],
                [60.4626, -34.1751, 39.4387],
                3.1819,
            ),
            (
                [63.0109, -31.0961, -5.8663],
                [62.8187, -29.7946, -4.0864],
                2.2133,
            ),
            (
                [61.2901, 3.7196, -5.3901],
                [61.4292, 2.2480, -4.9620],
                1.5389,
            ),
            (
                [35.0831, -44.1164, 3.7933],
                [35.0232, -40.0716, 1.5901],
                4.6063,
            ),
            (
                [22.7233, 20.0904, -46.6940],
                [23.0331, 14.9730, -42.5619],
                6.5847,
            ),
            (
                [36.4612, 47.8580, 18.3852],
                [36.2715, 50.5065, 21.2231],
                3.8864,
            ),
            (
                [90.8027, -2.0831, 1.4410],
                [91.1528, -1.6435, 0.0447],
                1.5051,
            ),
            (
                [90.9257, -0.5406, -0.9208],
                [88.6381, -0.8985, -0.7239],
                2.3238,
            ),
            (
                [6.7747, -0.2908, -2.4247],
                [5.8714, -0.0985, -2.2286],
                0.9441,
            ),
            (
                [2.0776, 0.0795, -1.1350],
                [0.9033, -0.0636, -0.5514],
                1.3191,
            ),
        ];

        let xyz_d65 = Cie1931.xyz_d65();
        for &(lab1_arr, lab2_arr, expected) in cases {
            let lab1 = CieLab::new(lab1_arr, xyz_d65);
            let lab2 = CieLab::new(lab2_arr, xyz_d65);
            let de = lab1.ciede(&lab2).unwrap();
            approx::assert_abs_diff_eq!(de, expected, epsilon = 1e-4);
        }
    }
}
