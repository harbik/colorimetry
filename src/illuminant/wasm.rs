// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

//! JS-WASM Interface code

#[cfg(feature = "cfi")]
use super::CFI;
#[cfg(feature = "cri")]
use super::CRI;
use super::{CieIlluminant, Illuminant};
use crate::spectrum::Spectrum;
use crate::traits::Light;
use wasm_bindgen::prelude::wasm_bindgen;

#[wasm_bindgen]
impl Illuminant {
    /// Create a new illuminant spectrum from the given data.
    ///
    /// The data must be the 401 values from 380 to 780 nm, with an interval size of 1 nanometer.
    #[wasm_bindgen(constructor)]
    pub fn new_js(data: &[f64]) -> Result<Illuminant, wasm_bindgen::JsError> {
        Ok(Illuminant(Spectrum::try_from(data)?))
    }

    /// Returns the spectral data values, as a Float64Array containing 401 data
    /// points, over a wavelength domain from 380 t0 780 nanometer, with a
    /// stepsize of 1 nanometer.
    #[wasm_bindgen(js_name=Values)]
    pub fn values_js(&self) -> Box<[f64]> {
        let values = self.spectrum().as_array().as_slice().to_vec();
        values.into_boxed_slice()
    }

    /// Calculates the Color Rendering Index for this illuminant spectrum.
    ///
    /// Returns a `CRI` object with a `ra()` method that gives the general colour
    /// rendering index Rₐ (average of R₁…R₈, scaled 0–100).
    #[cfg(feature = "cri")]
    #[wasm_bindgen(js_name=cri)]
    pub fn cri_js(&self) -> Result<crate::illuminant::CRI, crate::Error> {
        self.cri()
    }

    /// Calculates the Colour Fidelity Index (CFI / Rf) for this illuminant spectrum.
    ///
    /// Returns a `CFI` object exposing `rf()`, `rg()`, `rfHj()`, `rcsHj()`, `rhsHj()`,
    /// and `specialIndices()`.  Requires the `cfi` feature.
    #[cfg(feature = "cfi")]
    #[wasm_bindgen(js_name = cfi)]
    pub fn cfi_js(&self) -> Result<CFI, crate::Error> {
        self.cfi()
    }

    /// Get the CieIlluminant spectrum. Typically you don't need to use the Spectrum itself, as many
    /// methods just accept the CieIlluminant directly.
    #[wasm_bindgen(js_name=illuminant)]
    pub fn llluminant_js(stdill: CieIlluminant) -> Self {
        // need this as wasm_bindgen does not support `impl` on Enum types (yet?).
        // in Rust use CieIlluminant.spectrum() directly, which also gives a reference instead of a copy.
        stdill.illuminant().clone()
    }
}

/// Expose `CRI.ra()` to JavaScript when the `cri` feature is enabled.
#[cfg(feature = "cri")]
#[wasm_bindgen]
impl CRI {
    /// Returns the general colour rendering index Rₐ (0–100),
    /// the average of the first eight special rendering indices R₁…R₈.
    #[wasm_bindgen(js_name = ra)]
    pub fn ra_js(&self) -> f64 {
        self.ra()
    }
}

/// Expose `CFI` methods to JavaScript when the `cfi` feature is enabled.
#[cfg(feature = "cfi")]
#[wasm_bindgen]
impl CFI {
    /// General colour fidelity index Rf (0–100).
    ///
    /// A single overall score measuring how faithfully the test source renders the 99
    /// Colour Evaluation Samples compared to the reference illuminant.
    #[wasm_bindgen(js_name = colorFidelityIndex)]
    pub fn color_fidelity_index_js(&self) -> f64 {
        self.general_color_fidelity_index()
    }

    /// General colour gamut index Rg.
    ///
    /// Measures the area of the gamut polygon relative to the reference (100 = same area).
    /// Values above 100 indicate a wider gamut than the reference; below 100 means narrower.
    #[wasm_bindgen(js_name = colorGamutIndex)]
    pub fn color_gamut_index_js(&self) -> f64 {
        self.general_color_gamut_index()
    }

    /// Local colour fidelity index Rf,hj for each of the 16 hue bins (TM-30 / CIE 224:2017 §4.5).
    ///
    /// Returns a `Float64Array` of 16 values.  Bin 0 starts at 0° (red), progressing
    /// counter-clockwise in 22.5° steps around the hue circle.
    #[wasm_bindgen(js_name = localColorFidelityIndices)]
    pub fn local_color_fidelity_indices_js(&self) -> Box<[f64]> {
        Box::new(self.rf_hj())
    }

    /// Chroma shift index Rcs,hj for each of the 16 hue bins (TM-30 / CIE 224:2017 §4.6).
    ///
    /// Returns a `Float64Array` of 16 values.  Positive means the test source boosts
    /// saturation in that hue direction; negative means desaturation.  Typical range ≈ −0.5…+0.5.
    #[wasm_bindgen(js_name = chromaShiftIndices)]
    pub fn chroma_shift_indices_js(&self) -> Box<[f64]> {
        Box::new(self.rcs_hj())
    }

    /// Hue shift index Rhs,hj for each of the 16 hue bins (TM-30 / CIE 224:2017 §4.7).
    ///
    /// Returns a `Float64Array` of 16 values in radians, wrapped to (−π, π].
    /// Positive means a counter-clockwise hue shift; negative means clockwise.
    #[wasm_bindgen(js_name = hueShiftIndices)]
    pub fn hue_shift_indices_js(&self) -> Box<[f64]> {
        Box::new(self.rhs_hj())
    }

    /// Special colour fidelity indices Rf,i for all 99 CES (CIE 224:2017 §7).
    ///
    /// Returns a `Float64Array` of 99 values, one per Colour Evaluation Sample.
    #[wasm_bindgen(js_name = specialColorFidelityIndices)]
    pub fn special_color_fidelity_indices_js(&self) -> Box<[f64]> {
        Box::new(self.special_color_fidelity_indices())
    }
}
