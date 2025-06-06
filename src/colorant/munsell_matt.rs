//! Munsell Matt Spectral Data
//!
//! Data derived from data measured by a team from the University of Eastern Finland, using a Perkin-Elmer lambda 9 UV/VIS/NIR spectrofotometer,
//! for 1269 matt Munsell chips.
//! The original spectral data ranges from 380 to 800, with steps of 1nm, but here this dataset was reduced to a range from 380 to 780nm,
//! with steps of 5nm, by averaging, to improve calculation speed, and reduce program size.
//! These are measured data, for the specific Munsell chips, and by no means can be considered to be the nominal, or average spectral distribution
//! for all Munsell chips: please use it as an approximate representation, with unknown spectral deviation from the avarge reflection spectra of Munsell
//! chips.

mod data;

use std::{collections::BTreeMap, sync::LazyLock};

use crate::{
    cam::{CieCam16, ViewConditions},
    colorant::munsell_matt::data::{MUNSELL_MATT_DATA, MUNSELL_MATT_KEYS},
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
pub struct MunsellMatt(String, Spectrum);

impl MunsellMatt {
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
        let spectrum = Spectrum::linear_interpolate(
            &[
                *SPECTRUM_WAVELENGTH_RANGE.start() as f64,
                *SPECTRUM_WAVELENGTH_RANGE.end() as f64,
            ],
            &data,
        )
        .unwrap();
        MunsellMatt(key.to_string(), spectrum)
    }

    /// Try to find the spectral data for the given key.
    ///
    /// Fails if the key is not found.
    pub fn try_new(key: impl AsRef<str>) -> Result<MunsellMatt, Error> {
        if let Some(&i) = MM_KEY_MAP.get(key.as_ref()) {
            Ok(Self::new(i))
        } else {
            Err(Error::SpectrumNotFound(key.as_ref().to_string()))
        }
    }

    pub fn key(&self) -> &str {
        &self.0
    }

    pub fn spectrum(&self) -> &Spectrum {
        &self.1
    }
}

/// Collection of all Munsell Matt chips.
/// This collection can be iterated over to access each Munsell Matt chip's spectral data.
/// It provides a convenient way to access all available Munsell Matt chips for spectral analysis or color matching tasks.  
pub struct MunsellMattCollection;

impl MunsellMattCollection {
    /// Returns the number of Munsell Matt chips in the collection.
    pub fn len() -> usize {
        MATT_M
    }

    pub fn match_ciecam16(
        colorant: &dyn Filter,
        opt_illuminant: Option<&Illuminant>,
        opt_vc: Option<ViewConditions>,
        opt_observer: Option<observer::Observer>,
    ) -> Result<(String, f64), Error> {
        let illuminant = opt_illuminant.unwrap_or(&D65);
        let vc = opt_vc.unwrap_or_default();
        let observer = opt_observer.unwrap_or_default();

        let xyz = observer.xyz(illuminant, Some(colorant));
        let xyzn = observer.xyz(illuminant, None);
        let tgt_cam = CieCam16::from_xyz(xyz, xyzn, vc)?;
        let mut best_key = String::new();
        let mut best_delta_e = f64::MAX;
        for mm in MunsellMattCollection.into_iter() {
            let xyz_mm = observer.xyz(illuminant, Some(&mm));
            let cam_mm = CieCam16::from_xyz(xyz_mm, xyzn, vc)?;
            let delta_e = tgt_cam.de_ucs(&cam_mm)?;
            if delta_e < best_delta_e {
                best_delta_e = delta_e;
                best_key = mm.key().to_string();
            }
        }
        Ok((best_key, best_delta_e))
    }
}

pub struct MunsellMattIterator(usize);

impl Iterator for MunsellMattIterator {
    type Item = MunsellMatt;

    fn next(&mut self) -> Option<Self::Item> {
        if self.0 < MATT_M {
            let mm = MunsellMatt::new(self.0);
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

impl ExactSizeIterator for MunsellMattIterator {}

impl IntoIterator for MunsellMattCollection {
    type Item = MunsellMatt;

    type IntoIter = MunsellMattIterator;

    fn into_iter(self) -> Self::IntoIter {
        MunsellMattIterator(0)
    }
}

/// This implements the `Filter` trait, allowing it to be used in spectral calculations.
/// It provides access to the spectrum of a specific Munsell Matt chip, which can be used for color matching or other spectral analysis tasks.
impl Filter for MunsellMatt {
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
mod test_munsell_matt {
    use super::MunsellMattCollection;
    use crate::observer::Observer::Cie1931;
    #[test]
    fn test_iter() {
        MunsellMattCollection
            .into_iter()
            .enumerate()
            .for_each(|(i, m)| {
                let lab_d65 = Cie1931.lab_d65(&m);
                let key = m.0;
                println!("{i} {key} {lab_d65:?}");
            });

        let mm = crate::colorant::MunsellMattCollection
            .into_iter()
            .last()
            .unwrap();
        assert_eq!(mm.0, "10RP4/12".to_string());
    }

    #[test]
    fn test_match_ciecam16() {
        let colorant = crate::colorant::MunsellMatt::try_new("10RP4/12").unwrap();
        let (key, delta_e) = MunsellMattCollection::match_ciecam16(
            &colorant,
            None,
            None,
            Some(Cie1931),
        )
        .unwrap();
        assert_eq!(key, "10RP4/12");
        approx::assert_abs_diff_eq!(delta_e, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_match_r9() {
        let r9 = &crate::colorant::tcs::TCS[8];
        let (key, delta_e) = MunsellMattCollection::match_ciecam16(
            r9,
            None,
            None,
            None,
        )
        .unwrap();
        assert_eq!(key, "5R4/14");
        approx::assert_abs_diff_eq!(delta_e, 3.0, epsilon = 5e-2);
    }
}
