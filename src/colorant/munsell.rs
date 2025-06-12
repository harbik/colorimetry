//! Munsell Spectral Data
//!
//! This dataset contains spectral measurements of 1,269 matte Munsell color chips,
//! captured by a team from the University of Eastern Finland using a Perkin-Elmer Lambda 9 UV/VIS/NIR spectrophotometer.
//!
//! The original measurements span wavelengths from 380 nm to 800 nm at 1 nm intervals.
//! For improved performance and reduced data size, the dataset presented here has been downsampled to 5 nm intervals,
//! covering the range from 380 nm to 780 nm. This was done by averaging the original data.
//!
//! Please note that these are empirical measurements of specific physical Munsell chips,
//! and should not be considered nominal or representative of the entire Munsell color system.
//! Spectral deviations from the average reflectance of Munsell colors are unknown,
//! so use this dataset as an approximate reference only.

mod data;

use std::{collections::BTreeMap, sync::LazyLock};

use crate::{
    cam::{CieCam16, ViewConditions},
    colorant::munsell::data::{MUNSELL_MATT_DATA, MUNSELL_MATT_KEYS},
    error::Error,
    illuminant::D65,
    observer,
    prelude::Illuminant,
    spectrum::{Spectrum, SPECTRUM_WAVELENGTH_RANGE},
    traits::Filter,
};

pub(crate) const MATT_N: usize = 81;
pub(crate) const MATT_M: usize = 1269;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
/// A measured Munsell color sample, including its name and spectral reflectance data.
///
/// Each instance of `Munsell` represents one matte Munsell chip, identified by a string key
/// (e.g., `"5R 5/14"`), and its corresponding reflectance spectrum sampled at 1 nm intervals
/// from 380 nm to 780 nm. The spectral data is interpolated using Sprague’s method.
///
/// Note: These are empirical measurements of specific chips, and by no means represent typical
/// values. Please use this dataset as an approximate reference only.
///
pub struct Munsell(String, Spectrum);

impl Munsell {
    /// Get MunsellMatt Spectral data by index.
    ///
    /// Index value should be in the range from 0..1269. This will panic if a larger number is
    /// requested.
    pub fn new(i: usize) -> Self {
        let key = MUNSELL_MATT_KEYS[i];
        let data: [f64; MATT_N] = MUNSELL_MATT_DATA[i * MATT_N..(i + 1) * MATT_N]
            .iter()
            .map(|v| *v as f64)
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();
        let spectrum = Spectrum::sprague_interpolate(
            [
                *SPECTRUM_WAVELENGTH_RANGE.start() as f64,
                *SPECTRUM_WAVELENGTH_RANGE.end() as f64,
            ],
            &data,
        )
        .unwrap();
        Munsell(key.to_string(), spectrum)
    }

    /// Try to find the spectral data for the given key.
    ///
    /// Fails if the key is not found.
    pub fn try_new(key: impl AsRef<str>) -> Result<Munsell, Error> {
        if let Some(&i) = MM_KEY_MAP.get(key.as_ref()) {
            Ok(Self::new(i))
        } else {
            Err(Error::SpectrumNotFound(key.as_ref().to_string()))
        }
    }

    /// Returns the key for this Munsell chip.
    /// This key is a string identifier that uniquely represents the Munsell chip,
    /// typically in the format "Hue Value/Chroma" (e.g., "5R 5/14").
    pub fn key(&self) -> &str {
        &self.0
    }

    /// Returns the spectral reflectance data for this Munsell chip.
    /// The spectrum is represented as a `Spectrum` object, which contains the reflectance values
    /// at specific wavelengths (from 380 nm to 780 nm).
    /// This data can be used for color matching, spectral analysis, or other applications
    /// that require precise spectral information.
    pub fn spectrum(&self) -> &Spectrum {
        &self.1
    }
}

/// Collection of all Munsell Matt chips.
/// This collection can be iterated over to access each Munsell Matt chip's spectral data.
/// It provides a convenient way to access all available Munsell Matt chips for spectral analysis or color matching tasks.  
pub struct MunsellCollection;

impl MunsellCollection {
    /// Returns the number of Munsell Matt chips in the collection.
    pub fn len() -> usize {
        MATT_M
    }

    /// Finds the Munsell chip that best matches a given colorant using the CIECAM16-UCS color appearance model.
    ///
    /// This method compares the provided `colorant` against all matte Munsell samples and
    /// returns the key of the closest match along with the ΔE (color difference) value.
    ///
    /// The `colorant` must implement the [`Filter`] trait to provide spectral data.
    ///
    /// # Optional Parameters
    /// - `opt_illuminant`: The light source used for comparison. Defaults to D65 if `None`.
    /// - `opt_vc`: The viewing conditions (e.g., surround, luminance). Defaults to standard values if `None`.
    /// - `opt_observer`: The standard observer for XYZ conversion. Defaults to CIE 1931 if `None`.
    ///
    /// # Returns
    /// A `Result` containing:
    /// - `Ok((key, delta_e))`: the key of the closest Munsell chip and the corresponding color difference.
    /// - `Err`: if an error occurs during the color appearance transformation or ΔE computation.
    ///
    /// # Errors
    /// Returns an error if any step of the CIECAM16 transformation or UCS distance calculation fails.
    pub fn match_ciecam16(
        colorant: &dyn Filter,
        opt_illuminant: Option<&Illuminant>,
        opt_vc: Option<ViewConditions>,
        opt_observer: Option<observer::Observer>,
    ) -> Result<(String, f64), Error> {
        let illuminant = opt_illuminant.unwrap_or(&D65);
        let vc = opt_vc.unwrap_or_default();
        let observer = opt_observer.unwrap_or_default();
        let rxyz = observer.rel_xyz(illuminant, colorant);
        let tgt_cam = CieCam16::from_xyz(rxyz, vc);
        let mut best_key = String::new();
        let mut best_delta_e = f64::MAX;
        for mm in MunsellCollection.into_iter() {
            let xyz_mm = observer.rel_xyz(illuminant, &mm);
            let cam_mm = CieCam16::from_xyz(xyz_mm, vc);
            let delta_e = tgt_cam.de_ucs(&cam_mm)?;
            if delta_e < best_delta_e {
                best_delta_e = delta_e;
                best_key = mm.key().to_string();
            }
        }
        Ok((best_key, best_delta_e))
    }
}

pub struct MunsellIterator(usize);

impl Iterator for MunsellIterator {
    type Item = Munsell;

    fn next(&mut self) -> Option<Self::Item> {
        if self.0 < MATT_M {
            let mm = Munsell::new(self.0);
            self.0 += 1;
            Some(mm)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (MATT_M, Some(MATT_M))
    }
}

impl ExactSizeIterator for MunsellIterator {}

impl IntoIterator for MunsellCollection {
    type Item = Munsell;

    type IntoIter = MunsellIterator;

    fn into_iter(self) -> Self::IntoIter {
        MunsellIterator(0)
    }
}

/// This implements the `Filter` trait, allowing it to be used in spectral calculations.
/// It provides access to the spectrum of a specific Munsell Matt chip, which can be used for color matching or other spectral analysis tasks.
impl Filter for Munsell {
    fn spectrum(&self) -> std::borrow::Cow<'_, Spectrum> {
        (&self.1).into()
    }
}

/// A static map for Munsell Matt keys to their indices.
/// This map is used for quick lookups to find the index of a Munsell Matt by its key.
/// The keys are the Munsell Matt identifiers, and the values are their corresponding indices.
/// This map is immutable and can be used throughout the application for efficient access.
pub static MM_KEY_MAP: LazyLock<BTreeMap<&str, usize>> = LazyLock::new(|| {
    let mut map = BTreeMap::new();
    MUNSELL_MATT_KEYS
        .iter()
        .enumerate()
        .for_each(|(value, &key)| {
            map.insert(key, value);
        });
    map
});

#[cfg(test)]
mod test_munsell {
    use super::MunsellCollection;
    use crate::observer::Observer::Cie1931;
    #[test]
    fn test_iter() {
        MunsellCollection
            .into_iter()
            .enumerate()
            .for_each(|(i, m)| {
                let lab_d65 = Cie1931.lab_d65(&m);
                let key = m.0;
                println!("{i} {key} {lab_d65:?}");
            });

        let mm = crate::colorant::MunsellCollection
            .into_iter()
            .last()
            .unwrap();
        assert_eq!(mm.0, "10RP4/12".to_string());
    }

    #[test]
    fn test_match_ciecam16() {
        let colorant = crate::colorant::Munsell::try_new("10RP4/12").unwrap();
        let (key, delta_e) =
            MunsellCollection::match_ciecam16(&colorant, None, None, Some(Cie1931)).unwrap();
        assert_eq!(key, "10RP4/12");
        approx::assert_abs_diff_eq!(delta_e, 0.0, epsilon = 1e-6);
    }

    #[test]
    #[cfg(feature = "cri")]
    fn test_match_r9() {
        let r9 = &crate::colorant::tcs::TCS[8];
        let (key, delta_e) = MunsellCollection::match_ciecam16(r9, None, None, None).unwrap();
        assert_eq!(key, "5R4/14");
        approx::assert_abs_diff_eq!(delta_e, 2.8, epsilon = 5e-2);
    }
}
