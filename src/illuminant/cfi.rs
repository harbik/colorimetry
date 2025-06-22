use crate::math::distance;
use crate::observer::Observer::Cie1964;
use crate::xyz::RelXYZ;
use crate::{
    cam::{CamTransforms, CieCam02, TM30VC},
    colorant::{CES, N_CFI},
    illuminant::Illuminant,
};

use super::CCT;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
/// Container for CIE 2017 Colour Fidelity Index (**R<sub>f</sub>**) calculations,
/// including both the general color fidelity **R<sub>f</sub>** score and the 99 special color fidelity indices (**R<sub>f,1</sub>** to **R<sub>f,99</sub>**)
/// as specified in [CIE 224:2017](https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use).
///
/// # Requirements
/// - Requires the `cfi` feature to access color evaluation samples (CES) used for testing.
/// - Requires the `supplemental` feature to use the CIE 1964 observer, which is required in this model. Included when you enable the `cfi` feature.
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

    /// Returns the array of special color fidelity indices (Rf<sub>1</sub> through Rf<sub>99</sub>) as defined in
    /// [CIE 224:2017 – CIE 2017 Colour Fidelity Index for accurate scientific use](https://cie.co.at/publications/cie-2017-colour-fidelity-index-accurate-scientific-use).
    ///
    /// # Returns
    /// An array of 99 `f64` values. Each value represents the fidelity score (Rf,i) for the corresponding
    /// CES under the current light source compared to the reference illuminant (daylight or Planckian).
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
