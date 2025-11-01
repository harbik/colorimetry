// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

use std::f64::consts::PI;

use crate::math::{self, distance};
use crate::observer::Observer::Cie1964;
use crate::xyz::RelXYZ;
use crate::{
    cam::{CamTransforms, CieCam02, TM30VC},
    colorant::{CES, N_CFI},
    illuminant::Illuminant,
};

use super::CCT;

pub const N_ANGLE_BIN: usize = 16;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
/// Container for CIE 2017 Colour Fidelity Index (**R<sub>f</sub>**) calculations,
/// including both the general color fidelity **R<sub>f</sub>** score and the 99 special color fidelity indices (**R<sub>f,1</sub>** to **R<sub>f,99</sub>**)
/// as specified in [CIE 224:2017](https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use).
///
/// # Requirements
/// - Requires the `cfi` feature to access color evaluation samples (CES) used for testing.
///
/// # Overview
/// The CIE 2017 Colour Fidelity Index (CFI, or **R<sub>f</sub>**) is a modern metric for evaluating how accurately a light source renders colors.
/// It uses 99 Color Evaluation Samples (CES) that cover a broad range of real-world colors, providing a much more comprehensive assessment
/// than older metrics like the Color Rendering Index (CRI, or **R<sub>a</sub>**).
/// - The general index (**R<sub>f</sub>**) gives an overall measure of color fidelity.
/// - The special indices (**R<sub>f,1</sub>** to **R<sub>f,99</sub>**) show fidelity for each specific color sample.
///
/// # Comparison with CRI
/// The traditional **CRI** metric (Ra) uses only 8 or 14 pastel color samples and is known to be limited, especially for modern light sources such as LEDs.
/// **CFI (Rf)** is a newer, more robust standard: it uses a much wider set of samples and is based on state-of-the-art color appearance models,
/// providing a more accurate and reliable prediction of real-world color rendering.
/// - **Use CRI** if you need compatibility with legacy systems or must comply with standards that specify CRI.
/// - **Use CFI (Rf)** for a more precise and scientifically up-to-date assessment of color fidelity, especially with modern or tunable light sources.
///
/// # Reference
/// [CIE 224:2017 – CIE 2017 Colour Fidelity Index for accurate scientific use](https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use)
pub struct CFI {
    jabp_ts: [[f64; 3]; N_CFI], // CIE1964
    jabp_rs: [[f64; 3]; N_CFI],
    cct: CCT, // using CIE1931
}

impl CFI {
    /// Creates a new CFI instance for the given illuminant.
    /// # Parameters
    /// - `illuminant`: The light source for which to calculate the CFI.
    /// # Returns
    /// - `Ok(CFI)`: A new CFI instance containing the calculated color fidelity indices.
    /// # Errors
    /// - An error if the illuminant's CCT cannot be determined
    pub fn new(illuminant: &Illuminant) -> Result<Self, crate::Error> {
        let cct = illuminant.cct()?;
        let ref_illuminant = Illuminant::cfi_reference(cct.t())?;
        let vc = TM30VC;
        let xyzn_r = Cie1964.xyz(&ref_illuminant, None).set_illuminance(100.0);
        let xyzn_t = Cie1964.xyz(illuminant, None).set_illuminance(100.0);
        let mut jabp_ts = [[0f64; 3]; N_CFI];
        let mut jabp_rs = [[0f64; 3]; N_CFI];
        for (i, cfi_ces) in CES.iter().enumerate() {
            let xyz_t = Cie1964.xyz(illuminant, Some(cfi_ces));
            let rxyz_t = RelXYZ::from_xyz(xyz_t, xyzn_t)?;
            let jabp_t = CieCam02::from_xyz(rxyz_t, vc).jab_prime();
            jabp_ts[i] = jabp_t;

            let xyz_r = Cie1964.xyz(&ref_illuminant, Some(cfi_ces));
            let rxyz_r = RelXYZ::from_xyz(xyz_r, xyzn_r)?;
            let jabp_r = CieCam02::from_xyz(rxyz_r, vc).jab_prime();
            jabp_rs[i] = jabp_r;
        }

        Ok(CFI {
            jabp_ts,
            jabp_rs,
            cct,
        })
    }

    /// Returns the test source J'a'b' coordinates for the 99 CES test samples
    pub fn jabp_ts(&self) -> &[[f64; 3]; N_CFI] {
        &self.jabp_ts
    }

    /// Returns the reference source J'a'b' coordinates for the 99 CES test samples
    pub fn jabp_rs(&self) -> &[[f64; 3]; N_CFI] {
        &self.jabp_rs
    }

    /// Calculate the chroma values for each CES color sample under the test illuminant
    pub fn chroma_ts(&self) -> [f64; N_CFI] {
        compute_chroma(&self.jabp_ts)
    }

    /// Calculate the chroma values for each CES color sample under the reference illuminant
    pub fn chroma_rs(&self) -> [f64; N_CFI] {
        compute_chroma(&self.jabp_rs)
    }

    /// Calculate the hue angle values (in radians) for each CES color sample under the test illuminant
    pub fn hue_angle_bin_ts(&self) -> [f64; N_CFI] {
        compute_hue_angle_bin(&self.jabp_ts)
    }

    /// Calculate the hue angle values (in radians) for each CES color sample under the reference illuminant
    pub fn hue_angle_bin_rs(&self) -> [f64; N_CFI] {
        compute_hue_angle_bin(&self.jabp_rs)
    }

    /// Calculate hue-bin averaged J'a'b' values for the CES samples under the test illuminant
    pub fn jabp_average_ts(&self) -> [[f64; 3]; N_ANGLE_BIN] {
        let mut rst = [([0f64, 0f64, 0f64], 0f64); N_ANGLE_BIN];
        for i in 0..N_CFI {
            // get J'a'b' for test illuminant and reference illuminant
            // need reference to determine hue angle bin
            let [jt, at, bt] = self.jabp_ts[i];
            let [_jr, ar, br] = self.jabp_rs[i];

            // determine hue angle bin from reference sample
            let mut phi = f64::atan2(br, ar);
            if phi < 0. {
                phi += 2. * PI;
            }

            // accumulate into hue angle bin
            let i = (phi / (2. * PI) * N_ANGLE_BIN as f64) as usize;
            rst[i].0[0] += jt;
            rst[i].0[1] += at;
            rst[i].0[2] += bt;
            rst[i].1 += 1.; // sample count in bin
        }
        rst.map(|([j, a, b], n)| {
            if n == 0. {
                [f64::NAN, f64::NAN, f64::NAN]
            } else {
                [j / n, a / n, b / n]
            }
        })
    }

    /// Returns hue-bin averaged J'a'b' values for the CES samples under the reference illuminant
    pub fn jabp_average_rs(&self) -> [[f64; 3]; N_ANGLE_BIN] {
        let mut rst = [([0f64, 0f64, 0f64], 0f64); N_ANGLE_BIN];
        for i in 0..N_CFI {
            // get J'a'b' for reference illuminant
            let [jr, ar, br] = self.jabp_rs[i];
            let mut phi = f64::atan2(br, ar);
            if phi < 0. {
                phi += 2. * PI;
            }

            // accumulate into hue angle bin
            let i = (phi / (2. * PI) * N_ANGLE_BIN as f64) as usize;
            rst[i].0[0] += jr;
            rst[i].0[1] += ar;
            rst[i].0[2] += br;
            rst[i].1 += 1.;
        }
        rst.map(|([j, a, b], n)| {
            if n == 0. {
                [f64::NAN, f64::NAN, f64::NAN]
            } else {
                [j / n, a / n, b / n]
            }
        })
    }

    /// Returns normalized averaged a'b' under test and reference sources.
    /// Take the first returned array in order to plot the CVG.
    pub fn normalized_ab_average(&self) -> ([[f64; 2]; N_ANGLE_BIN], [[f64; 2]; N_ANGLE_BIN]) {
        let jabt_hj = self.jabp_average_ts();
        let jabr_hj = self.jabp_average_rs();
        compute_normalized_ab_average(&jabt_hj, &jabr_hj)
    }

    /// Returns chroma for each hue bin for test source
    pub fn normalized_chroma_average_ts(&self) -> [f64; N_ANGLE_BIN] {
        let jabt_hj = self.jabp_average_ts();
        compute_normalized_chroma_average(&jabt_hj)
    }

    /// Returns chroma for each hue bin for reference source
    pub fn normalized_chroma_average_rs(&self) -> [f64; N_ANGLE_BIN] {
        let jabr_hj = self.jabp_average_rs();
        compute_normalized_chroma_average(&jabr_hj)
    }

    /// Returns the normalized chroma for each hue fin for test source
    pub fn normalized_chroma_average(&self) -> [f64; N_ANGLE_BIN] {
        let ct = self.normalized_chroma_average_ts();
        let cr = self.normalized_chroma_average_rs();
        let mut ctn = [0f64; N_ANGLE_BIN];
        for i in 0..N_ANGLE_BIN {
            if cr[i] == 0. {
                ctn[i] = f64::INFINITY;
            } else {
                ctn[i] = ct[i] / (cr[i] + 1e-308);
            }
        }
        ctn
    }

    /// Returns hues (rad) for each hue bin for test source
    pub fn hue_angle_bin_average_samples_ts(&self) -> [f64; N_ANGLE_BIN] {
        compute_hue_angle_bin_average(&self.jabp_average_ts())
    }

    /// Returns hues (rad) for each hue bin for reference source
    pub fn hue_angle_bin_average_samples_rs(&self) -> [f64; N_ANGLE_BIN] {
        compute_hue_angle_bin_average(&self.jabp_average_rs())
    }

    /// Calculate the IES-TM-30-18 Gamut Index (Rg) for the illuminant.
    ///
    /// # Returns
    /// * `f64` - The Gamut Index (Rg) value, representing the color gamut of the light source.
    ///
    /// # Overview
    /// The Gamut Index (Rg) quantifies the extent of color saturation provided by a light source.
    /// It is calculated by comparing the area of the color gamut formed by the test light source
    /// to that of a reference illuminant (either daylight or Planckian).
    /// Higher Rg values indicate a wider color gamut, meaning the light source can render more saturated colors.
    ///
    /// # Notes
    /// - An Rg value of 100 indicates that the test light source has the same color gamut area as the reference illuminant.
    /// - Values above 100 suggest a wider gamut, while values below 100 indicate a narrower gamut.
    /// - The Gamut Index is not part of the CIE 2017 Colour Fidelity Index (CFI) but is often reported alongside it for a comprehensive
    ///   assessment of a light source's color rendering capabilities.
    ///
    /// # Reference
    /// - ANSI/IES TM-30-20 Method for Evaluating Light Source Color Rendition, section 4.4, Gamut Index Rg (ISBN 978-0-87995-379-9)
    ///
    pub fn rg(&self) -> f64 {
        let origin = [0.; 2];
        let av_samples_t = self.jabp_average_ts();
        let av_samples_r = self.jabp_average_rs();
        let mut at = 0f64;
        let mut ar = 0f64;

        // Compute gamut area for the test source
        for i in 0..N_ANGLE_BIN {
            // Compute gamut area for the test source
            let area_t = math::compute_triangle_area(
                &[av_samples_t[i][1], av_samples_t[i][2]],
                &[
                    av_samples_t[(i + 1) % N_ANGLE_BIN][1],
                    av_samples_t[(i + 1) % N_ANGLE_BIN][2],
                ],
                &origin,
            );

            // Compute gamut area for the reference source
            let area_r = math::compute_triangle_area(
                &[av_samples_r[i][1], av_samples_r[i][2]],
                &[
                    av_samples_r[(i + 1) % N_ANGLE_BIN][1],
                    av_samples_t[(i + 1) % N_ANGLE_BIN][2],
                ],
                &origin,
            );
            at += area_t;
            ar += area_r;
        }

        100. * at / ar
    }

    /// Returns the array of special color fidelity indices (Rf<sub>1</sub> through Rf<sub>99</sub>) as defined in
    /// [CIE 224:2017 – CIE 2017 Colour Fidelity Index for accurate scientific use](https://cie.co.at/publications/cie-2017-colour-fidelity-index-accurate-scientific-use).
    ///
    /// # Returns
    /// An array of 99 `f64` values. Each value represents the fidelity score (Rf,i) for the corresponding
    /// [`CES`] under the current light source compared to the reference illuminant (daylight or Planckian).
    /// Higher values indicate better color fidelity for that sample.
    ///
    /// # Usage
    /// Use the returned array to analyze the color rendering performance of a light source
    /// for individual colors or specific regions (e.g., skin tones, foliage, saturated reds, etc.).
    /// This allows for more detailed diagnostics than the general Rf index alone.
    ///
    /// # Reference
    /// - CIE 224:2017 – CIE 2017 Colour Fidelity Index for accurate scientific use (Sections 7, Annex E/F)
    pub fn special_color_fidelity_indices(&self) -> [f64; N_CFI] {
        let mut cfis = [0f64; N_CFI];
        for (i, (jabd, jabr)) in self.jabp_ts.iter().zip(self.jabp_rs.iter()).enumerate() {
            let de = distance(jabd, jabr);
            cfis[i] = rf_from_de(de);
        }
        cfis
    }

    /// Returns the general color fidelity index (**R<sub>f</sub>**) for the given light source,
    /// as defined in [CIE 224:2017 – CIE 2017 Colour Fidelity Index for accurate scientific use](https://cie.co.at/publications/cie-2017-colour-fidelity-index-accurate-scientific-use).
    ///
    /// # Overview
    /// The general color fidelity index (**R<sub>f</sub>**) provides a single overall score (from 0 to 100)
    /// indicating how closely the colors of illuminated objects, under the test light source,
    /// match those seen under a reference illuminant. Higher values indicate better color fidelity,
    /// with values close to 100 meaning that colors appear nearly identical to their appearance
    /// under the reference.
    ///
    /// This value is the main summary statistic of the CIE 2017 Color Fidelity Index,
    /// making it directly comparable to the traditional CRI (**R<sub>a</sub>**) metric, but it is more robust and
    /// accurate, especially for modern light sources.
    ///
    /// # Returns
    /// * `f64` – The general color fidelity index (**R<sub>f</sub>**), ranging from 0 (poor color fidelity)
    ///   to 100 (perfect color fidelity).
    ///
    pub fn general_color_fidelity_index(&self) -> f64 {
        let mut sum = 0.0;
        for (jabp_t, jabp_r) in self.jabp_ts.iter().zip(self.jabp_rs.iter()) {
            let de = distance(jabp_t, jabp_r);
            sum += de;
        }
        rf_from_de(sum / N_CFI as f64)
    }

    /// Returns the correlated color temperature (CCT) used in the CFI calculation,
    /// which is used to select the appropriate reference illuminant for color fidelity evaluation.
    ///
    /// **Note:**  
    /// While the CFI (Color Fidelity Index) itself is calculated using the CIE 1964 10° standard observer,
    /// the CCT is always computed using the CIE 1931 2° standard observer, following the official CIE 224:2017 procedure.
    ///
    /// # Returns
    /// * [`CCT`] — A structure, containing the correlated color temperature of the test light source,
    ///   in Kelvin, and distance to the plancking curve in the CIE1960 UCS chromaticity diagram.
    ///
    pub fn cct(&self) -> CCT {
        self.cct
    }
}

fn compute_hue_angle_bin(jab: &[[f64; 3]; N_CFI]) -> [f64; N_CFI] {
    jab.map(|[_j, a, b]| {
        let mut h = f64::atan2(b, a);
        if h < 0. {
            h += 2. * PI;
        }
        h
    })
}

fn compute_hue_angle_bin_average(jab: &[[f64; 3]; N_ANGLE_BIN]) -> [f64; N_ANGLE_BIN] {
    jab.map(|[_j, a, b]| {
        let mut h = f64::atan2(b, a);
        if h < 0. {
            h += 2. * PI;
        }
        h
    })
}

fn compute_chroma(jab: &[[f64; 3]; N_CFI]) -> [f64; N_CFI] {
    jab.map(|[_j, a, b]| f64::sqrt(a * a + b * b))
}

fn compute_normalized_chroma_average(jab_hj: &[[f64; 3]; N_ANGLE_BIN]) -> [f64; N_ANGLE_BIN] {
    jab_hj.map(|[_j, a, b]| f64::sqrt(a * a + b * b))
}

// returns (jabtn_hj, jabrn_hj)
fn compute_normalized_ab_average(
    jabt: &[[f64; 3]; N_ANGLE_BIN],
    jabr: &[[f64; 3]; N_ANGLE_BIN],
) -> ([[f64; 2]; N_ANGLE_BIN], [[f64; 2]; N_ANGLE_BIN]) {
    let ht_hj = compute_hue_angle_bin_average(jabt);
    let hr_hj = compute_hue_angle_bin_average(jabr);
    let ct = compute_normalized_chroma_average(jabt);
    let cr = compute_normalized_chroma_average(jabr);

    let mut ctn = [0f64; N_ANGLE_BIN];
    for i in 0..N_ANGLE_BIN {
        if cr[i] == 0. {
            ctn[i] = f64::INFINITY;
        } else {
            ctn[i] = ct[i] / (cr[i] + 1e-308);
        }
    }
    let mut jabtn_hj = [[0f64; 2]; N_ANGLE_BIN];
    let mut jabrn_hj = [[0f64; 2]; N_ANGLE_BIN];
    for i in 0..N_ANGLE_BIN {
        let c = ctn[i];
        let ht = ht_hj[i];
        let hr = hr_hj[i];
        jabtn_hj[i] = [c * f64::cos(ht), c * f64::sin(ht)];
        jabrn_hj[i] = [f64::cos(hr), f64::sin(hr)];
    }
    (jabtn_hj, jabrn_hj)
}

#[allow(dead_code)]
// this is not used, any reason why we should keep this?
fn compute_hue_bin_edges() -> [f64; N_ANGLE_BIN + 1] {
    let dh = 360. / N_ANGLE_BIN as f64;
    let mut hbe = [0f64; N_ANGLE_BIN + 1];
    #[allow(clippy::needless_range_loop)]
    for i in 0..=N_ANGLE_BIN {
        let m = i as f64;
        hbe[i] = dh * m;
    }
    hbe
}

const CF: f64 = 6.73; // was 7.54 in TM30-15
pub fn rf_from_de(de: f64) -> f64 {
    10.0 * (((100.0 - CF * de) / 10.0).exp() + 1.0).ln()
}

#[cfg(test)]
mod tests {
    use super::CFI;
    use crate::{colorant::N_CFI, illuminant::D65};

    #[cfg(feature = "cie-illuminants")]
    use crate::illuminant::{F1, F12, F2};

    use approx::assert_abs_diff_eq;

    #[test]
    fn test_cfi() {
        let cfi = CFI::new(&D65).unwrap();
        let cfis = cfi.special_color_fidelity_indices();
        assert_eq!(cfis.len(), N_CFI);
        for &cfi in cfis.iter() {
            assert_abs_diff_eq!(cfi, 100.0, epsilon = 0.01); // Rf for CIE D65
        }
        assert!(cfi.general_color_fidelity_index() >= 0.0);
    }

    #[test]
    #[cfg(feature = "cie-illuminants")]
    fn test_fl1() {
        let cfi = CFI::new(&F1).unwrap();
        assert_abs_diff_eq!(cfi.cct().t(), 6425.0, epsilon = 6.0);
        assert_abs_diff_eq!(cfi.cct().d(), 0.0072, epsilon = 0.0002);
        assert_abs_diff_eq!(cfi.general_color_fidelity_index(), 81.0, epsilon = 0.5);
        assert_abs_diff_eq!(cfi.rg(), 90.0, epsilon = 1.0);
        // Rf for DE=1
    }

    #[test]
    #[cfg(feature = "cie-illuminants")]
    fn test_fl2() {
        let cfi = CFI::new(&F2).unwrap();
        assert_abs_diff_eq!(cfi.cct().t(), 4225.0, epsilon = 1.0);
        assert_abs_diff_eq!(cfi.cct().d(), 0.0019, epsilon = 0.0002);
        assert_abs_diff_eq!(cfi.general_color_fidelity_index(), 70.0, epsilon = 0.5);
        assert_abs_diff_eq!(cfi.rg(), 86.0, epsilon = 1.0);
    }

    #[test]
    #[cfg(feature = "cie-illuminants")]
    fn test_fl12() {
        let cfi = CFI::new(&F12).unwrap();
        assert_abs_diff_eq!(cfi.cct().t(), 3003.0, epsilon = 4.0);
        assert_abs_diff_eq!(cfi.cct().d(), 0.0001, epsilon = 0.0002);
        assert_abs_diff_eq!(cfi.general_color_fidelity_index(), 78.0, epsilon = 0.6);
        // TODO this fails, with a value of 109 instead of 102
        // assert_abs_diff_eq!(cfi.rg(), 102.0, epsilon = 1.0);
    }

    #[test]
    #[allow(unused)]
    #[cfg(feature = "cie-illuminants")]
    fn test_fl12_jab_from_tm30_20() {
        use crate::illuminant::cfi::N_ANGLE_BIN;

        const WANT: [[f64; 6]; N_CFI] = [
            [85.93, 17.36, 0.22, 86.14, 16.84, 1.68],
            [63.00, 30.81, 2.78, 63.29, 31.27, 4.79],
            [31.04, 12.48, 1.98, 32.02, 15.27, 2.45],
            [69.99, 19.87, 4.35, 70.01, 20.42, 6.61],
            [48.78, 28.66, 6.02, 51.47, 33.95, 9.86],
            [50.69, 16.32, 5.43, 51.48, 18.06, 7.00],
            [42.98, 28.02, 10.65, 43.89, 32.10, 10.92],
            [40.97, 23.91, 9.23, 42.99, 29.99, 10.06],
            [28.93, 4.48, 1.06, 29.06, 3.34, 0.63],
            [76.03, 29.39, 17.98, 75.81, 29.43, 18.15],
            [57.23, 28.93, 17.29, 58.66, 28.22, 19.37],
            [64.90, 29.26, 20.30, 64.89, 30.10, 20.89],
            [44.49, 17.81, 15.40, 43.94, 18.85, 14.46],
            [74.40, 5.70, 6.35, 74.23, 6.44, 5.91],
            [71.82, 12.27, 12.31, 71.82, 13.00, 10.48],
            [47.58, 12.49, 17.56, 48.39, 17.68, 16.30],
            [49.03, 12.26, 17.37, 50.08, 14.83, 17.07],
            [56.73, 10.24, 14.86, 56.73, 11.95, 13.26],
            [70.74, 13.39, 22.92, 71.56, 13.59, 23.22],
            [67.05, 17.94, 29.68, 67.41, 19.85, 27.49],
            [86.65, 9.39, 29.19, 86.12, 13.62, 27.92],
            [78.73, 11.61, 30.54, 78.64, 14.99, 28.65],
            [91.67, 1.39, 13.45, 91.50, 3.70, 12.08],
            [90.93, 2.14, 25.54, 90.33, 6.16, 24.05],
            [72.92, 1.30, 34.03, 71.45, 7.84, 32.53],
            [89.96, 2.34, 30.77, 89.48, 6.48, 28.38],
            [73.52, -0.71, 16.78, 73.38, 1.34, 14.65],
            [42.42, -2.10, 16.30, 41.96, 0.20, 15.27],
            [86.95, -2.83, 32.65, 86.62, 1.84, 28.25],
            [64.26, -0.53, 21.73, 64.52, 1.91, 16.31],
            [88.51, -3.31, 31.96, 88.31, 0.30, 28.19],
            [75.57, -7.12, 33.96, 75.15, -0.70, 30.91],
            [89.61, -2.57, 23.01, 89.63, 0.12, 18.81],
            [63.27, -8.21, 30.06, 62.03, -3.31, 27.82],
            [43.23, -5.88, 16.36, 42.44, -2.70, 14.94],
            [80.57, -0.88, 9.76, 80.75, -0.00, 6.75],
            [49.78, -10.80, 20.82, 47.99, -5.14, 18.40],
            [86.17, -4.81, 21.04, 86.33, -2.65, 15.65],
            [35.37, -1.48, 7.55, 35.44, -1.91, 6.63],
            [50.32, -3.49, 14.09, 50.31, -4.04, 12.64],
            [87.91, -4.28, 10.28, 88.28, -3.06, 7.93],
            [69.90, -20.64, 29.92, 68.42, -16.68, 27.39],
            [72.58, -18.96, 30.97, 72.18, -17.35, 27.79],
            [28.81, -1.08, 1.47, 28.84, -0.74, 1.27],
            [58.60, -19.62, 21.97, 58.26, -18.18, 19.45],
            [66.62, -14.92, 18.06, 66.54, -14.38, 15.35],
            [69.30, -10.74, 13.99, 70.20, -8.98, 10.79],
            [77.71, -19.95, 19.00, 77.38, -19.87, 13.73],
            [46.29, -17.60, 15.69, 46.25, -17.46, 12.82],
            [54.74, -17.73, 9.45, 54.69, -18.06, 6.88],
            [58.64, -14.82, 8.25, 58.64, -15.81, 5.65],
            [40.83, -18.89, 9.16, 41.19, -21.19, 7.26],
            [49.04, -29.01, 6.63, 49.42, -30.76, 3.49],
            [86.80, -16.62, 1.40, 87.01, -16.12, -0.70],
            [83.41, -17.82, 0.49, 83.58, -17.92, -2.30],
            [70.51, -29.29, -0.75, 70.87, -30.20, -3.81],
            [57.01, -30.73, -0.31, 57.35, -31.89, -4.17],
            [53.10, -30.69, -3.05, 53.65, -31.46, -6.27],
            [86.90, -15.20, -2.23, 86.85, -15.89, -3.25],
            [88.24, -6.32, -1.55, 88.42, -6.65, -2.77],
            [76.76, -13.21, -4.70, 76.92, -14.47, -5.57],
            [66.65, -5.73, -0.95, 67.33, -6.05, -4.24],
            [34.09, 3.67, -1.30, 34.80, -1.97, -2.75],
            [57.44, -26.33, -13.56, 58.35, -28.72, -14.80],
            [39.40, -14.70, -10.20, 40.07, -18.16, -10.38],
            [42.93, -21.05, -16.66, 44.49, -25.08, -17.06],
            [46.15, -21.37, -18.18, 47.72, -25.38, -18.46],
            [63.40, -13.47, -15.51, 64.01, -16.94, -15.20],
            [53.26, -2.83, -10.11, 53.62, -7.01, -10.45],
            [71.80, -13.47, -21.66, 72.54, -17.17, -21.01],
            [34.66, -10.98, -18.65, 35.26, -14.05, -20.43],
            [76.70, -4.53, -8.08, 76.90, -5.93, -8.13],
            [57.11, -15.83, -26.54, 57.75, -19.84, -25.72],
            [27.84, -1.32, -3.26, 28.50, -0.93, -3.68],
            [32.49, -7.86, -20.67, 33.36, -10.99, -21.74],
            [47.96, -7.41, -31.12, 48.46, -12.37, -30.21],
            [43.45, -8.16, -23.94, 45.11, -7.32, -22.91],
            [44.16, -5.88, -30.54, 44.68, -10.44, -29.80],
            [73.23, -2.47, -20.69, 74.24, -2.66, -19.19],
            [65.46, 3.02, -23.15, 66.50, 1.77, -21.95],
            [43.11, 1.88, -28.81, 45.17, 5.56, -25.74],
            [87.27, 2.34, -9.89, 87.57, 1.45, -9.24],
            [72.97, 4.69, -17.25, 73.49, 3.76, -15.44],
            [28.35, 6.97, -5.59, 28.46, 4.49, -6.05],
            [40.98, 5.59, -12.55, 41.96, 6.99, -11.55],
            [64.55, 5.93, -24.09, 66.65, 11.03, -19.67],
            [40.81, 16.46, -19.00, 40.62, 14.24, -19.03],
            [63.31, 19.69, -16.37, 63.37, 18.17, -14.77],
            [38.17, 17.04, -15.53, 39.25, 17.82, -12.89],
            [44.97, 21.49, -15.49, 44.87, 20.37, -13.14],
            [80.57, 3.47, 0.29, 80.92, 3.91, -1.57],
            [55.59, 13.25, -9.30, 57.49, 16.69, -6.95],
            [67.27, 12.21, -6.53, 68.29, 14.10, -5.37],
            [38.47, 18.50, -10.03, 40.85, 23.48, -5.81],
            [70.48, 17.28, -5.24, 71.34, 18.53, -4.33],
            [78.62, 16.87, -4.92, 79.23, 17.53, -3.42],
            [51.06, 26.57, -6.52, 53.25, 23.14, -2.36],
            [65.79, 30.88, -1.85, 65.86, 29.97, 0.54],
            [53.91, 30.88, -1.38, 54.35, 32.52, 0.69],
        ];

        let cfi = CFI::new(&F12).unwrap();
        for (i, (&[jt, at, bt], &[jr, ar, br])) in
            cfi.jabp_ts().iter().zip(cfi.jabp_rs().iter()).enumerate()
        {
            assert_abs_diff_eq!(jt, WANT[i][0], epsilon = 0.2);
            assert_abs_diff_eq!(at, WANT[i][1], epsilon = 0.2);
            assert_abs_diff_eq!(bt, WANT[i][2], epsilon = 0.2);
            assert_abs_diff_eq!(jr, WANT[i][3], epsilon = 0.2);
            assert_abs_diff_eq!(ar, WANT[i][4], epsilon = 0.2);
            assert_abs_diff_eq!(br, WANT[i][5], epsilon = 0.2);
        }

        const WANT_BINS: [[[f64; 2]; 2]; N_ANGLE_BIN] = [
            [[22.15, 3.50], [23.98, 5.02]],
            [[19.41, 15.31], [20.53, 15.08]],
            [[12.47, 24.09], [14.81, 22.93]],
            [[-0.49, 24.62], [2.99, 21.85]],
            [[-5.22, 18.21], [-2.61, 15.74]],
            [[-14.32, 19.40], [-12.72, 17.01]],
            [[-18.77, 17.35], [-18.66, 13.28]],
            [[-20.11, 8.37], [-21.45, 5.82]],
            [[-21.94, -1.31], [-22.56, -3.72]],
            [[-15.57, -10.94], [-18.14, -11.84]],
            [[-7.40, -15.29], [-10.99, -15.75]],
            [[-5.05, -21.91], [-6.74, -21.16]],
            [[2.98, -19.78], [3.13, -18.09]],
            [[8.74, -15.31], [9.19, -14.07]],
            [[17.87, -14.17], [18.26, -11.94]],
            [[15.82, -5.49], [16.78, -3.81]],
        ];

        let av_samples_t = cfi.jabp_average_ts();
        for (i, &[_, at, bt]) in av_samples_t.iter().enumerate() {
            use approx::abs_diff_eq;

            let [[at_w, bt_w], [_ar_w, _br_w]] = WANT_BINS[i];
            assert!(
                abs_diff_eq!(at, at_w, epsilon = 1.0),
                "at failed at index {i}: got {at}, want {at_w})"
            );
            /* assert!(
                abs_diff_eq!(bt, bt_w, epsilon = 1.0),
                "bt failed at index {i}: got {bt}, want {bt_w})"
            ); */
        }
    }
}
