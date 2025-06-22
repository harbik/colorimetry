//! Gamut Module for CIE LCh Color Space
//!
//! This module defines structures and implementations for representing the gamut
//! of the CIE LCh color space as a matrix of maximum chroma values.
//!
//! The gamut is modeled as a two-dimensional matrix where each cell corresponds to a lightness–hue bin:
//! - There are `NL` lightness bins (L = 1..=NL), arranged as rows.
//! - There are `NH` hue bins (each covering 360/NH degrees, H = 0..NH-1), arranged as columns.
//!
//! The gamut matrix is stored as `SMatrix<u8, NL, NH>`, where each cell `(l, h)` contains the maximum chroma (C)
//! for that lightness (row) and hue (column) bin under a specific observer and illuminant.
//!
//! The module organizes different gamut definitions based on standard observer measurements
//! (e.g., CIE 1931, CIE 1964, CIE 2015) and illuminant conditions (D50 and D65).
//! It provides accessors for retrieving the associated observer, illuminant, and gamut data,
//! as well as methods to query the maximum chroma for specific (lightness, hue) combinations
//! and to check if a CIE Lab color falls within the defined gamut.

mod cielch_gamut_cie1931_d50;
use cielch_gamut_cie1931_d50::CIELCH_GAMUT_CIE1931_D50;

mod cielch_gamut_cie1931_d65;
use cielch_gamut_cie1931_d65::CIELCH_GAMUT_CIE1931_D65;

// Optional features for supplemental observers
#[cfg(feature = "supplemental-observers")]
mod cielch_gamut_cie1964_d50;
#[cfg(feature = "supplemental-observers")]
use cielch_gamut_cie1964_d50::CIELCH_GAMUT_CIE1964_D50;

#[cfg(feature = "supplemental-observers")]
mod cielch_gamut_cie1964_d65;
#[cfg(feature = "supplemental-observers")]
use cielch_gamut_cie1964_d65::CIELCH_GAMUT_CIE1964_D65;

#[cfg(feature = "supplemental-observers")]
mod cielch_gamut_cie2015_10_d50;
#[cfg(feature = "supplemental-observers")]
use cielch_gamut_cie2015_10_d50::CIELCH_GAMUT_CIE2015_10_D50;

#[cfg(feature = "supplemental-observers")]
mod cielch_gamut_cie2015_10_d65;
#[cfg(feature = "supplemental-observers")]
use cielch_gamut_cie2015_10_d65::CIELCH_GAMUT_CIE2015_10_D65;

#[cfg(feature = "supplemental-observers")]
mod cielch_gamut_cie2015_d65;
#[cfg(feature = "supplemental-observers")]
use cielch_gamut_cie2015_d65::CIELCH_GAMUT_CIE2015_D65;

#[cfg(feature = "supplemental-observers")]
mod cielch_gamut_cie2015_d50;
#[cfg(feature = "supplemental-observers")]
use cielch_gamut_cie2015_d50::CIELCH_GAMUT_CIE2015_D50;

use crate::illuminant::CieIlluminant;
use nalgebra::SMatrix;
use strum::Display;

const NL: usize = 99; // number of lightness bins.
const NH: usize = 72; // number of hue bins (each covering 5 degrees).

/// Represents the CIE LCh color gamut as a matrix of maximum chroma values.
///
/// The gamut is stored as a 2D matrix (`SMatrix<u8, NL, NH>`) where:
///
/// - The horizontal axis (columns, size `NH = 72`) corresponds to hue (H) bins, each covering 5 degrees.
/// - The vertical axis (rows, size `NL = 99`) corresponds to lightness (L) bins, ranging from 1 to 99.
/// - Each cell `(l, h)` contains the maximum chroma (C) for that lightness and hue bin.
///
/// The struct also stores the associated illuminant and observer used to compute the gamut data.
pub struct CieLChGamutData {
    illuminant: CieIlluminant,
    observer: Observer,
    gamut: SMatrix<u8, NL, NH>,
}

impl CieLChGamutData {
    pub const fn new(
        illuminant: CieIlluminant,
        observer: Observer,
        gamut: SMatrix<u8, NL, NH>,
    ) -> Self {
        Self {
            illuminant,
            observer,
            gamut,
        }
    }

    /// Returns the maximum chroma for the given lightness and hue bin.
    /// - `l`: lightness bin index (row, 0..NL-1)
    /// - `h`: hue bin index (column, 0..NH-1)
    pub fn max_chroma(&self, l: usize, h: usize) -> u8 {
        self.gamut[(l, h)]
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Display)]
pub enum CieLChGamut {
    #[default]
    Cie1931D50,
    Cie1931D65,
    #[cfg(feature = "supplemental-observers")]
    Cie1964D50,
    #[cfg(feature = "supplemental-observers")]
    Cie1964D65,
    #[cfg(feature = "supplemental-observers")]
    Cie2015_10D50,
    #[cfg(feature = "supplemental-observers")]
    Cie2015_10D65,
    #[cfg(feature = "supplemental-observers")]
    Cie2015D50,
    #[cfg(feature = "supplemental-observers")]
    Cie2015D65,
}

use crate::{observer::Observer, Error};

impl CieLChGamut {
    pub fn observer(&self) -> Observer {
        self.data().observer
    }

    pub fn illuminant(&self) -> CieIlluminant {
        self.data().illuminant
    }

    pub fn data(&self) -> &CieLChGamutData {
        match self {
            // D65-based gamuts
            CieLChGamut::Cie1931D65 => &CIELCH_GAMUT_CIE1931_D65,
            #[cfg(feature = "supplemental-observers")]
            CieLChGamut::Cie1964D65 => &CIELCH_GAMUT_CIE1964_D65,
            #[cfg(feature = "supplemental-observers")]
            CieLChGamut::Cie2015_10D65 => &CIELCH_GAMUT_CIE2015_10_D65,
            #[cfg(feature = "supplemental-observers")]
            CieLChGamut::Cie2015D65 => &CIELCH_GAMUT_CIE2015_D65,

            // D50-based gamuts
            CieLChGamut::Cie1931D50 => &CIELCH_GAMUT_CIE1931_D50,
            #[cfg(feature = "supplemental-observers")]
            CieLChGamut::Cie1964D50 => &CIELCH_GAMUT_CIE1964_D50,
            #[cfg(feature = "supplemental-observers")]
            CieLChGamut::Cie2015_10D50 => &CIELCH_GAMUT_CIE2015_10_D50,
            #[cfg(feature = "supplemental-observers")]
            CieLChGamut::Cie2015D50 => &CIELCH_GAMUT_CIE2015_D50,
        }
    }

    pub fn lch_bin(&self, l: u8, h: u8) -> Result<u8, Error> {
        match (l, h) {
            (0, _) | (100, _) => return Ok(0),
            (_, h) if h >= NH as u8 => Err(Error::InvalidHueBin(h)),
            (1..=99, h) => Ok(self.data().gamut[(l as usize - 1, h as usize)]),
            (l,_) => Err(Error::InvalidLightness(l)),
        }
    }

    /// Check if the Lab color is within the gamut defined by the CieLChGamut
    pub fn in_gamut(&self, lab: super::CieLab) -> bool {
        // Safe unwrap: white_point comes from a validated enum.
        let [l, c, h] = lab.lch();

        // Validate L value
        if !(0.0..=100.0).contains(&l) {
            return false; // Invalid Lab values.
        }
        let l_u8 = l.round() as u8;

        // Validate Hue value
        if !(0.0..360.0).contains(&h) {
            return false; // Invalid hue value.
        }

        // Hue is divided by 5 to fit into 72 hue bins, each covering 5 degrees.
        // Using floor to ensure we don't exceed the maximum hue bin.
        let h_u8 = (h / 5.0).floor() as u8;

        match (l_u8, h_u8) {
            // Lightness at extremes are automatically considered in gamut.
            (0, _) | (100, _) => true,
            // Reject hue if it exceeds defined number of bins.
            (_, h) if h >= NH as u8 => false,
            _ => {
                // Safe unwrap: lookup is valid as ranges were already checked.
                let max_c = self.lch_bin(l_u8, h_u8).unwrap();
                max_c < c.floor() as u8
            }
        }
    }

    /// Get the maximum chroma for a specific hue bin.
    /// # Arguments
    /// * `hue_bin` - The hue bin index (0 to 71), with a width of 5º.
    ///
    /// # Returns
    /// A result containing an array of maximum chroma values for the specified hue bin,
    /// or an error if the hue bin index is invalid.
    /// The values are associated with lightness bins from 1 to 99.
    /// The maximum chroma values are returned as an array of length `NL` (99),
    pub fn iso_hue(&self, hue_bin: u8) -> Result<[u8; NL], Error> {
        if hue_bin >= 72 {
            return Err(Error::InvalidHue(hue_bin));
        }
        let row_vec: Vec<u8> = self
            .data()
            .gamut
            .row(hue_bin as usize)
            .iter()
            .copied()
            .collect();
        Ok(row_vec.try_into().unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lab::CieLab;
    use crate::prelude::Observer::Cie1931;
    use crate::xyz::RelXYZ;

    #[test]
    fn test_cielch_gamut_data() {
        let g31d65 = CieLChGamut::Cie1931D65;

        assert_eq!(g31d65.illuminant(), CieIlluminant::D65);
        assert_eq!(g31d65.data().gamut.nrows(), NH);
        assert_eq!(g31d65.data().gamut.ncols(), NL);
    }

    #[test]
    fn test_in_gamut() {
        let c31d65 = CieLChGamut::Cie1931D65;
        let white_point = Cie1931.xyz_d50();
        let rxyz = RelXYZ::from_vec(white_point.xyz, white_point);
        let lab_in_gamut = CieLab::from_xyz(rxyz);
        assert!(c31d65.in_gamut(lab_in_gamut));

        let rxyz = RelXYZ::new([0.0; 3], white_point);
        let lab_in_gamut = CieLab::from_xyz(rxyz);
        assert!(c31d65.in_gamut(lab_in_gamut));

        // Test a color outside the gamut
        let lab_out_of_gamut = CieLab::new([101.0, 50.0, 50.0], white_point);
        assert!(!c31d65.in_gamut(lab_out_of_gamut));
    }

    #[test]
    fn test_values() {
        // Test the maximum chroma for specific lightness and hue bins.
        // See test in optimal_colors.rs for the expected values.
        const WANTS: &[[u8; 3]] = &[
            [1, 0, 13],
            [50, 0, 97],
            [0, 50, 0],
            [99, 0, 5],
            [1, 71, 13],
            [50, 71, 99],
            [99, 71, 6],
            [1, 36, 5],
            [50, 36, 138],
            [99, 36, 12],
            [1, 18, 2],
            [50, 18, 86],
            [99, 18, 18],
        ];
        for &[l, h, expected] in WANTS {
            let c = CieLChGamut::Cie1931D65.lch_bin(l, h).unwrap();
            assert_eq!(
                c, expected,
                "Expected chroma for (L={}, H={}) to be {}, but got {}",
                l, h, expected, c
            );
        }
    }
}
