/*!
# Color Rendering Index Calculation

# References
- CIE 013.3-1995 Method of measuring and specifying colour rendering properties of light sources

 */

use crate::{
    colorant::{N_TCS, TCS},
    error::Error,
    illuminant::Illuminant,
    observer::CIE1931,
    xyz::XYZ,
};

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Debug, Clone, Copy)]
/// The **Color Rendering Index (CRI)** for a light source, computed according to CIE 13.3-1995.
///
/// This struct holds the 14 individual rendering indices R₁…R₁₄ for the standard test color samples,
/// and provides the general CRI, Rₐ, which is the average of the first eight Rᵢ values.  
///
/// # Calculation Method
/// 1. The test illuminant is scaled to 100 lx and converted to CIE XYZ under the CIE 1931 observer.  
/// 2. Each of the 14 standard Colorant test spectra (TCS) is measured under both the test and the  
///    reference illuminant (black-body or D-series at the test’s correlated color temperature).  
/// 3. For each sample, the color difference ΔE in CIE UVW space is computed, and  
///    Rᵢ = 100 − 4.6 · ΔE.  
/// 4. The general CRI Rₐ is then  
///    ```text
///    Rₐ = (R₁ + R₂ + … + R₈) / 8
///    ```
///
/// # Examples
/// ```rust
/// use colorimetry::illuminant::{Illuminant, CRI};
///
/// // Compute CRI for the D65 illuminant:
/// let cri: CRI = (&Illuminant::d65()).try_into().unwrap();
///
/// // General CRI:
/// let ra = cri.ra();
/// println!("General CRI Rₐ = {:.1}", ra);
///
/// // All 14 individual Rᵢ values:
/// let ri_values = cri.values();
/// for (i, &ri) in ri_values.iter().enumerate() {
///     println!("R{} = {:.1}", i + 1, ri);
/// }
/// ```
///
/// # Notes
/// - This implementation uses the **CIE 1931** color space and requires the `"cri"` feature to be enabled in the crate.  
/// - The CRI-metric is now considered somewhat outdated; newer metrics (e.g., TM-30) are recommended for modern lighting applications.  
///   However, CRI remains widely used and understood across the lighting industry.
///
/// # Errors
/// Constructing a `CRI` can fail if the illuminant’s correlated color temperature is out of the
/// valid range (1000–25000 K) or its distance from the Planckian locus exceeds 0.05 Δuv.
///
pub struct CRI([f64; N_TCS]);

impl CRI {
    pub fn new(s: impl AsRef<Illuminant>) -> Result<Self, Error> {
        s.as_ref().try_into()
    }

    pub fn ra(&self) -> f64 {
        self.0.iter().take(8).sum::<f64>() / 8.0
    }

    pub fn values(&self) -> [f64; N_TCS] {
        self.0
    }
}

impl std::ops::Index<usize> for CRI {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

/// CRI calculation.
///
/// Can fail, for example if the Spectrum's correlated color temperature is out of range.
/// Uses CIE1931, and requires "cct"-feature.
impl TryFrom<&Illuminant> for CRI {
    type Error = Error;

    fn try_from(illuminant: &Illuminant) -> Result<Self, Self::Error> {
        let illuminant = &illuminant.clone().set_illuminance(&CIE1931, 100.0);
        // Calculate Device Under Test (dut) XYZ illuminant and sample values
        let xyz_dut = CIE1931.xyz_from_spectrum(illuminant.as_ref());
        let xyz_dut_samples: [XYZ; N_TCS] = TCS
            .each_ref()
            .map(|colorant| CIE1931.xyz(illuminant, Some(colorant)));

        // Determine reference color temperarture value
        let cct_dut = xyz_dut.cct()?.t();
        //println!("cct dut {cct_dut}");
        let illuminant_ref = if cct_dut <= 5000.0 {
            Illuminant::planckian(cct_dut).set_illuminance(&CIE1931, 100.0)
        } else {
            Illuminant::d_illuminant(cct_dut)?.set_illuminance(&CIE1931, 100.0)
        };

        // Calculate the reference illuminant values
        let xyz_ref = CIE1931.xyz_from_spectrum(illuminant_ref.as_ref());
        let xyz_ref_samples: [XYZ; N_TCS] = TCS
            .each_ref()
            .map(|colorant| CIE1931.xyz(&illuminant_ref, Some(colorant)));

        let cdt = cd(xyz_dut.uv60());
        let cdr = cd(xyz_ref.uv60());

        let ri: [f64; N_TCS] = xyz_ref_samples
            .iter()
            .zip(xyz_dut_samples.iter())
            .map(|(xyzr, xyz)| {
                let cdti = cd(xyz.uv60());
                let uv_vk = uv_kries(cdt, cdr, cdti);
                let xyz_vk = XYZ::from_luv60(uv_vk[0], uv_vk[1], Some(xyz.xyz.y), None).unwrap();
                let uvw = xyz_vk.uvw64(xyz_ref);
                let uvwr = xyzr.uvw64(xyz_ref);
                100.0
                    - 4.6
                        * ((uvw[0] - uvwr[0]).powi(2)
                            + (uvw[1] - uvwr[1]).powi(2)
                            + (uvw[2] - uvwr[2]).powi(2))
                        .sqrt()
            })
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();

        Ok(CRI(ri))
    }
}

impl AsRef<[f64]> for CRI {
    fn as_ref(&self) -> &[f64] {
        &self.0
    }
}

fn cd(uv60: [f64; 2]) -> [f64; 2] {
    let [u, v] = uv60;
    [
        (4.0 - u - 10.0 * v) / v,
        (1.708 * v - 1.481 * u + 0.404) / v,
    ]
}

fn uv_kries(cdt: [f64; 2], cdr: [f64; 2], cdti: [f64; 2]) -> [f64; 2] {
    let [ct, dt] = cdt;
    let [cr, dr] = cdr;
    let [cti, dti] = cdti;
    let den = 16.518 + 1.481 * (cr / ct) * cti - (dr / dt) * dti;
    [
        (10.872 + 0.404 * (cr / ct) * cti - 4.0 * (dr / dt) * dti) / den,
        5.520 / den,
    ]
}

#[cfg(test)]
mod cri_test {
    use crate::illuminant::{CRI, D50};

    #[test]
    fn cri_d50() {
        // should be all 100.0
        let cri0: CRI = (&D50).try_into().unwrap();
        // println!("{cri0:?}");
        approx::assert_ulps_eq!(
            cri0.as_ref(),
            [100.0; crate::illuminant::cri::N_TCS].as_ref(),
            epsilon = 0.03
        );
    }

    #[test]
    #[cfg(feature = "cie-illuminants")]
    fn cri_f1() {
        use crate::illuminant::CieIlluminant;
        // should be all 100.0

        let cri0: CRI = CieIlluminant::F1.illuminant().try_into().unwrap();
        println!("{cri0:?}");
        // approx::assert_ulps_eq!(
        //     cri0.as_ref(),
        //     [100.0; crate::illuminant::cri::N_TCS].as_ref(),
        //     epsilon = 0.05
        // );
    }

    #[test]
    #[cfg(feature = "cie-illuminants")]
    fn cri_f3_1() {
        use crate::illuminant::CieIlluminant;
        // 2932K, check with values as given in CIE15:2004 Table T.8.2
        let cri0: CRI = CieIlluminant::F3_1.illuminant().try_into().unwrap();
        approx::assert_ulps_eq!(
            cri0.as_ref(),
            [42, 69, 89, 39, 41, 52, 66, 13, -109, 29, 19, 21, 47, 93]
                .map(|v| v as f64)
                .as_ref(),
            epsilon = 1.0
        );
    }

    #[test]
    #[cfg(feature = "cie-illuminants")]
    fn cri_f3_11() {
        use crate::illuminant::CieIlluminant;
        // 5854K, check with values as given in CIE15:2004 Table T.8.2
        let cri0: CRI = CRI::new(CieIlluminant::F3_11).unwrap();

        let expected_ra = 78.0;
        approx::assert_abs_diff_eq!(cri0.ra(), expected_ra, epsilon = 1.0);

        let expected_values = [
            90.0, 86.0, 49.0, 82.0, 81.0, 70.0, 85.0, 79.0, 24.0, 34.0, 64.0, 50.0, 90.0, 67.0,
        ];

        approx::assert_abs_diff_eq!(
            cri0.values().as_ref(),
            expected_values.as_ref(),
            epsilon = 1.0
        );
    }
}
