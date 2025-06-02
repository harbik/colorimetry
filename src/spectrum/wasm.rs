//! JS-WASM Interface code

use super::Spectrum;
use wasm_bindgen::prelude::wasm_bindgen;

#[wasm_bindgen]
impl Spectrum {
    /// Create a new spectrum from the given data.
    ///
    /// The data must be the 401 values from 380 to 780 nm, with an interval size of 1 nanometer.
    ///
    /// If the Spectral data you have uses another wavelength domain and/or a different
    /// wavelength interval, use the linear interpolate constructor,
    /// which takes a wavelength domain and spectral data as arguments.
    #[wasm_bindgen(constructor)]
    pub fn new_js(data: &[f64]) -> Result<Spectrum, wasm_bindgen::JsError> {
        Ok(Spectrum::try_from(data)?)
    }

    /// Returns the spectral data values, as a Float64Array containing 401 data
    /// points, over a wavelength domain from 380 t0 780 nanometer, with a
    /// stepsize of 1 nanometer.
    #[wasm_bindgen(js_name=Values)]
    pub fn values_js(&self) -> Box<[f64]> {
        let values: &[f64] = self.as_ref();
        values.into()
    }

    /// This function maps spectral data with irregular intervals or intervals different than 1
    /// nanometer to the standard spectrum as used in this library.
    ///
    /// For domains with a regular interval, the wavelength slice should have a size of two, containing
    /// the minimum and maximum wavelength values, both also in units of meters or nanometers.
    ///
    /// For irregular domains, this function requires a slice of wavelengths and a slice of spectral
    /// data, both of the same size. The wavelengths can be specified in units of meters or nanometers.
    ///
    /// In case of duplicate wavelength values the last data values is used, so it is impossible to
    /// define filters with vertical edges using this method.
    ///
    /// ```rust
    /// // Creates a linear gradient filter, with a zero transmission at 380 nanometer, and full
    /// // transmission at 780 nanometer. This is an example using a uniform wavelength domain as input.
    /// use colorimetry::prelude::*;
    /// use approx::assert_ulps_eq;
    /// let data = [0.0, 1.0];
    /// let wl = [380.0, 780.0];
    /// let mut spd = Spectrum::linear_interpolate(&wl, &data).unwrap();
    /// assert_ulps_eq!(spd[380], 0.);
    /// assert_ulps_eq!(spd[380+100], 0.25);
    /// assert_ulps_eq!(spd[380+200], 0.5);
    /// assert_ulps_eq!(spd[380+300], 0.75);
    /// assert_ulps_eq!(spd[380+400], 1.0);
    ///
    /// // Creates a top hat filter, with slanted angles, using an irregular
    /// // wavelength domain.
    /// let data = vec![0.0, 1.0, 1.0, 0.0];
    /// let wl = vec![480.0, 490.0, 570.0, 580.0];
    /// let spd = Spectrum::linear_interpolate(&wl, &data).unwrap();
    /// assert_ulps_eq!(spd[380+0], 0.0);
    /// assert_ulps_eq!(spd[380+100], 0.0);
    /// assert_ulps_eq!(spd[380+110], 1.0);
    /// assert_ulps_eq!(spd[380+190], 1.0);
    /// assert_ulps_eq!(spd[380+200], 0.0);
    /// assert_ulps_eq!(spd[380+300], 0.0);
    /// assert_ulps_eq!(spd[380+400], 0.0);
    /// ```
    #[wasm_bindgen(js_name=linearInterpolate)]
    pub fn linear_interpolate_js(
        wavelengths: &[f64],
        data: &[f64],
    ) -> Result<Spectrum, crate::Error> {
        Self::linear_interpolate(wavelengths, data)
    }
}
