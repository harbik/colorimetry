// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

//! JS-WASM Interface code

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

    /// Calculates the Color Rendering Index values for illuminant spectrum.
    #[cfg(feature = "cri")]
    #[wasm_bindgen(js_name=cri)]
    pub fn cri_js(&self) -> Result<crate::illuminant::CRI, crate::Error> {
        self.cri()
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
