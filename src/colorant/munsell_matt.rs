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

use js_sys::Iter;
use std::{collections::BTreeMap, sync::LazyLock};
use wasm_bindgen::prelude::*;

use crate::{
    colorant::munsell_matt::data::{MUNSELL_MATT_DATA, MUNSELL_MATT_KEYS},
    error::CmtError,
    spectrum::{Spectrum, SPECTRUM_WAVELENGTH_RANGE},
    traits::Filter,
};

pub(crate) const MATT_N: usize = 81;
pub(crate) const MATT_M: usize = 1269;

#[wasm_bindgen]
pub struct MunsellMatt(String, Spectrum);

impl MunsellMatt {
    /// Get MunsellMatt Spectral data b`y index.
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
    pub fn try_new(key: String) -> Result<MunsellMatt, CmtError> {
        if let Some(&i) = MM_KEY_MAP.get(key.as_str()) {
            Ok(Self::new(i))
        } else {
            Err(CmtError::SpectrumNotFound(key))
        }
    }
}

// JS-WASM Interface code
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl MunsellMatt {}

pub struct MunsellMattCollection;

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

impl Filter for MunsellMatt {
    fn spectrum(&self) -> std::borrow::Cow<'_, Spectrum> {
        (&self.1).into()
    }
}

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
    use crate::prelude::*;

    #[test]
    fn test_iter() {
        /*
        MunsellMattCollection.into_iter().enumerate().for_each(|(i, m)| {
            let lab_d65 = CIE1931.lab_d65(&m);
            let key = m.0;
            println!("{i} {key} {lab_d65:?}");
        });
         */
        let mm = MunsellMattCollection.into_iter().last().unwrap();
        assert_eq!(mm.0, "10RP4/12".to_string());
    }
}
