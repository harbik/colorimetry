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

    /// Returns the test source J'a'b' coordinates for the 99 CES
    pub fn jabp_ts(&self) -> &[[f64; 3]; N_CFI] {
        &self.jabp_ts
    }

    /// Returns the reference source J'a'b' coordinates for the 99 CES
    pub fn jabp_rs(&self) -> &[[f64; 3]; N_CFI] {
        &self.jabp_rs
    }

    /// Returns chroma of each sample under the test source
    pub fn chroma_ts(&self) -> [f64; N_CFI] {
        compute_chroma(&self.jabp_ts)
    }

    /// Returns chroma of each sample under the refenrence source
    pub fn chorma_rs(&self) -> [f64; N_CFI] {
        compute_chroma(&self.jabp_rs)
    }

    /// Returns hue angles (rad) for each sample under the test source
    pub fn hue_angle_bin_ts(&self) -> [f64; N_CFI] {
        compute_hue_angle_bin(&self.jabp_ts)
    }

    /// Returns hue angles (rad) for each sample under the reference source
    pub fn hue_angle_bin_rs(&self) -> [f64; N_CFI] {
        compute_hue_angle_bin(&self.jabp_rs)
    }

    /// Returns hue-bin averaged J'a'b' under the test source
    pub fn jabp_average_ts(&self) -> [[f64; 3]; N_ANGLE_BIN] {
        let mut rst = [([0f64, 0f64, 0f64], 0f64); N_ANGLE_BIN];
        for i in 0..N_CFI {
            let [jt, at, bt] = self.jabp_ts[i];
            let [_jr, ar, br] = self.jabp_rs[i];
            let mut phi = f64::atan2(br, ar);
            if phi < 0. {
                phi = 2. * PI + phi;
            }
            let i = (phi / (2. * PI) * N_ANGLE_BIN as f64) as usize;
            rst[i].0[0] += jt;
            rst[i].0[1] += at;
            rst[i].0[2] += bt;
            rst[i].1 += 1.;
        }
        rst.map(|([j, a, b], n)| 
            if n == 0. {[f64::NAN, f64::NAN, f64::NAN]}
            else {[j / n, a / n, b / n]}
        )
    }

    /// Returns hue-bin averaged J'a'b' under the reference source
    pub fn jabp_average_rs(&self) -> [[f64; 3]; N_ANGLE_BIN] {
        let mut rst = [([0f64, 0f64, 0f64], 0f64); N_ANGLE_BIN];
        for i in 0..N_CFI {
            let [jr, ar, br] = self.jabp_rs[i];
            let mut phi = f64::atan2(br, ar);
            if phi < 0. {
                phi = 2. * PI + phi;
            }
            let i = (phi / (2. * PI) * N_ANGLE_BIN as f64) as usize;
            rst[i].0[0] += jr;
            rst[i].0[1] += ar;
            rst[i].0[2] += br;
            rst[i].1 += 1.;
        }
        rst.map(|([j, a, b], n)| 
            if n == 0. {[f64::NAN, f64::NAN, f64::NAN]}
            else {[j / n, a / n, b / n]}
        )
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

    /// Returns the gamut index
    pub fn rg(&self) -> f64 {
        let origin = [0.; 2];
        let av_samples_t = self.jabp_average_ts();
        let av_samples_r = self.jabp_average_rs();
        let mut at = 0f64;
        let mut ar = 0f64;
        // Compute area of At
        for i in 0..N_ANGLE_BIN {
            let area_t = math::compute_triangle_area(
                &[av_samples_t[i][1], av_samples_t[i][2]],
                &[av_samples_t[(i + 1) % N_ANGLE_BIN][1], av_samples_t[(i + 1) % N_ANGLE_BIN][2]],
                &origin
            );
            let area_r = math::compute_triangle_area(
                &[av_samples_r[i][1], av_samples_r[i][2]],
                &[av_samples_r[(i + 1) % N_ANGLE_BIN][1], av_samples_t[(i + 1) % N_ANGLE_BIN][2]],
                &origin
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
    /// * [`CCT`] — A strcuture, containing the correlated color temperature of the test light source,
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
            h = 2. * PI + h;
        } 
        h
    })
}

fn compute_hue_angle_bin_average(jab: &[[f64; 3]; N_ANGLE_BIN]) -> [f64; N_ANGLE_BIN] {
    jab.map(|[_j, a, b]| {
        let mut h = f64::atan2(b, a);
        if h < 0. {
            h = 2. * PI + h;
        } 
        h
    })
}

fn compute_chroma(jab: &[[f64; 3]; N_CFI]) -> [f64; N_CFI] {
    jab.map(|[_j, a, b]| {
        f64::sqrt(a * a + b * b)
    })
}

fn compute_normalized_chroma_average(jab_hj: &[[f64; 3]; N_ANGLE_BIN]) -> [f64; N_ANGLE_BIN] {
    jab_hj.map(|[_j, a, b]| {
        f64::sqrt(a * a + b * b)
    })
}

// returns (jabtn_hj, jabrn_hj)
fn compute_normalized_ab_average(jabt: &[[f64; 3]; N_ANGLE_BIN], jabr: &[[f64; 3]; N_ANGLE_BIN]) -> ([[f64; 2]; N_ANGLE_BIN], [[f64; 2]; N_ANGLE_BIN]) {
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

fn compute_hue_bin_edges() -> [f64; N_ANGLE_BIN + 1] {
    let dh = 360. / N_ANGLE_BIN as f64;
    let mut hbe = [0f64; N_ANGLE_BIN + 1];
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
        assert_abs_diff_eq!(cfi.cct().t(), 6428.0, epsilon = 3.0); // Rf for DE=1
        assert_abs_diff_eq!(cfi.general_color_fidelity_index(), 81.0, epsilon = 0.5);
        // Rf for DE=1
    }

    #[test]
    #[cfg(feature = "cie-illuminants")]
    fn test_fl2() {
        let cfi = CFI::new(&F2).unwrap();
        assert_abs_diff_eq!(cfi.cct().t(), 4225.0, epsilon = 1.0); // Rf for DE=1
        assert_abs_diff_eq!(cfi.general_color_fidelity_index(), 70.0, epsilon = 0.5);
    }

    #[test]
    #[cfg(feature = "cie-illuminants")]
    fn test_fl12() {
        let cfi = CFI::new(&F12).unwrap();
        assert_abs_diff_eq!(cfi.cct().t(), 3000.0, epsilon = 1.0); // Rf for DE=1
        assert_abs_diff_eq!(cfi.general_color_fidelity_index(), 78.0, epsilon = 0.6);
        // Rf for DE=1
    }
}
