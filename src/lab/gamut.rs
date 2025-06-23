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
use approx::abs_diff_eq;
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
/// - The vertical axis (rows, size `NL = 99`) corresponds to lightness (L) bins, ranging from 1 to 99.
/// - The horizontal axis (columns, size `NH = 72`) corresponds to hue (H) bins, each covering 5 degrees.
/// - Each cell `(l, h)` contains the maximum chroma (C) for that lightness (row) and hue (column) bin.
///
/// The struct also stores the associated illuminant and observer used to compute the gamut data.
pub struct CieLChGamutData {
    illuminant: CieIlluminant,
    observer: Observer,
    max_chroma_matrix: SMatrix<u8, NL, NH>,
}

impl CieLChGamutData {
    pub const fn new(
        illuminant: CieIlluminant,
        observer: Observer,
        max_chroma_matrix: SMatrix<u8, NL, NH>,
    ) -> Self {
        Self {
            illuminant,
            observer,
            max_chroma_matrix,
        }
    }

    /// Returns the maximum chroma for the given lightness and hue bin.
    /// - `l`: lightness bin index (row, 0..NL-1, corresponding to L = 1..=NL)
    /// - `h`: hue bin index (column, 0..NH-1, each covering 5 degrees)
    pub fn max_chroma(&self, l: usize, h: usize) -> u8 {
        self.max_chroma_matrix[(l, h)]
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

    pub fn max_chroma_bin(&self, l: u8, h: u8) -> Result<u8, Error> {
        match (l, h) {
            (0, _) | (100, _) => Ok(0),
            (_, h) if h >= NH as u8 => Err(Error::InvalidHueBin(h)),
            (1..=99, h) => Ok(self.data().max_chroma_matrix[(l as usize - 1, h as usize)]),
            (l, _) => Err(Error::InvalidLightnessBin(l)),
        }
    }

    /// Get the maximum chroma for a specific lightness and hue,
    /// by interpolation of the gamut data.
    /// # Arguments
    /// * `lightness` - The lightness value (0.0 to 100.0).
    /// * `hue` - The hue value (0.0 to 360.0).
    /// # Returns
    /// A result containing the maximum chroma value for the specified lightness and hue,
    /// or an error if the lightness or hue is out of bounds.
    /// # Errors
    /// Returns an error if the lightness or hue is out of bounds.
    /// The lightness must be in the range [0.0, 100.0] and the hue must be in the range [0.0, 360.0].
    pub fn max_chroma(&self, lightness: f64, hue: f64) -> Result<f64, crate::Error> {
        if abs_diff_eq!(lightness, 0.0, epsilon = 1e-6)
            || abs_diff_eq!(lightness, 100.0, epsilon = 1e-6)
        {
            return Ok(0.0); // Chroma is zero at lightness extremes.
        }
        if !(0.0..=100.0).contains(&lightness) {
            return Err(Error::InvalidLightness(lightness));
        }
        if !(0.0..360.0).contains(&hue) {
            return Err(Error::InvalidHue(hue));
        }
        // Convert lightness to bin indices
        let l_bin_low = lightness.floor() as u8;
        let l_bin_high = l_bin_low + 1;
        let h_bin_low = (hue / 5.0).floor() as u8;
        let h_bin_high = h_bin_low + 1;

        // Fractions for interpolation
        let l_frac = lightness - l_bin_low as f64;
        let h_frac = (hue / 5.0) - h_bin_low as f64;

        // Get the four surrounding chroma values
        let c00 = self.max_chroma_bin(l_bin_low, h_bin_low)? as f64;
        let c01 = self.max_chroma_bin(l_bin_low, h_bin_high)? as f64;
        let c10 = self.max_chroma_bin(l_bin_high, h_bin_low)? as f64;
        let c11 = self.max_chroma_bin(l_bin_high, h_bin_high)? as f64;

        // Bilinear interpolation
        let c0 = c00 * (1.0 - h_frac) + c01 * h_frac;
        let c1 = c10 * (1.0 - h_frac) + c11 * h_frac;
        let c = c0 * (1.0 - l_frac) + c1 * l_frac;

        Ok(c)
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
                let max_c = self.max_chroma_bin(l_u8, h_u8).unwrap();
                (c.floor() as u8) <= max_c
            }
        }
    }

    /// Get the maximum chroma for a specific hue bin (column).
    ///
    /// # Arguments
    /// * `hue_bin` - The hue bin index (column, 0 to 71), each covering 5º.
    ///
    /// # Returns
    /// An array of maximum chroma values for the specified hue bin,
    /// indexed by lightness bin (row, 0 to 98, corresponding to L = 1..=99).
    pub fn iso_hue(&self, hue_bin: u8) -> Result<[u8; NL], Error> {
        if hue_bin >= NH as u8 {
            return Err(Error::InvalidHueBin(hue_bin));
        }
        let mut chromas = [0; NL];
        let matrix = &self.data().max_chroma_matrix;
        for l in 0..NL {
            chromas[l] = matrix[(l, hue_bin as usize)];
        }
        Ok(chromas)
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
        assert_eq!(g31d65.data().max_chroma_matrix.nrows(), NL);
        assert_eq!(g31d65.data().max_chroma_matrix.ncols(), NH);
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
            let c = CieLChGamut::Cie1931D65.max_chroma_bin(l, h).unwrap();
            assert_eq!(
                c, expected,
                "Expected chroma for (L={}, H={}) to be {}, but got {}",
                l, h, expected, c
            );
        }
    }

    #[test]
    fn test_iso_hue() {
        let c31d65 = CieLChGamut::Cie1931D65;
        let hue_bin = 3; // Example hue bin
        let max_chroma = c31d65.iso_hue(hue_bin).unwrap();

        let want: Vec<u8> = CieLChGamut::Cie1931D65
            .data()
            .max_chroma_matrix
            .column(hue_bin as usize)
            .iter()
            .copied()
            .collect();

        // Check if the returned chroma values match the expected values
        for (i, &c) in max_chroma.iter().enumerate() {
            assert_eq!(
                c, want[i],
                "Mismatch at index {}: expected {}, got {}",
                i, want[i], c
            );
        }
    }

    #[test]
    fn test_max_chroma_at_bin_centers() {
        let gamut = CieLChGamut::Cie1931D65;
        // At bin center (L=50, H=100), should match the bin value exactly
        let l = 50.0;
        let h = 100.0;
        let expected = gamut.max_chroma_bin(l as u8, (h / 5.0) as u8).unwrap() as f64;
        let result = gamut.max_chroma(l, h).unwrap();
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_max_chroma_interpolation() {
        let gamut = CieLChGamut::Cie1931D65;
        // Between bins: halfway between L=50 and L=51, H=100 and H=105
        let l = 50.5;
        let h = 102.5;
        let c00 = gamut.max_chroma_bin(50, 20).unwrap() as f64;
        let c01 = gamut.max_chroma_bin(50, 21).unwrap() as f64;
        let c10 = gamut.max_chroma_bin(51, 20).unwrap() as f64;
        let c11 = gamut.max_chroma_bin(51, 21).unwrap() as f64;
        let l_frac = 0.5;
        let h_frac = 0.5;
        let c0 = c00 * (1.0 - h_frac) + c01 * h_frac;
        let c1 = c10 * (1.0 - h_frac) + c11 * h_frac;
        let expected = c0 * (1.0 - l_frac) + c1 * l_frac;
        let result = gamut.max_chroma(l, h).unwrap();
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_max_chroma_out_of_bounds() {
        let gamut = CieLChGamut::Cie1931D65;
        assert!(gamut.max_chroma(-1.0, 0.0).is_err());
        assert!(gamut.max_chroma(1.0, -1.0).is_err());
        assert!(gamut.max_chroma(101.0, 0.0).is_err());
        assert!(gamut.max_chroma(1.0, 360.0).is_err());
    }

    #[test]
    fn test_extrema() {
        let gamut = CieLChGamut::Cie1931D65;
        // Lightness at 0 should return 0 chroma
        assert_eq!(gamut.max_chroma(0.0, 180.0).unwrap(), 0.0);
        // Lightness at 100 should return 0 chroma
        assert_eq!(gamut.max_chroma(100.0, 180.0).unwrap(), 0.0);
    }

    #[test]
    fn test_low_lightness() {
        let gamut = CieLChGamut::Cie1931D65;
        // Lightness at 0.5, and hue 0.0, should return 6.5 chroma
        // Interpolate between 0 and first value in the first column of the gamut table
        approx::assert_abs_diff_eq!(gamut.max_chroma(0.5, 0.0).unwrap(), 6.5, epsilon = 1e-6);
        // Lightness at 0.5, and hue 0.5*5 = 2.5, should return 6.5 chroma as well
        // Interpolate between 0 and first value in the first column of the gamut table
        approx::assert_abs_diff_eq!(gamut.max_chroma(0.5, 2.5).unwrap(), 6.5, epsilon = 1e-6);
    }

    #[test]
    fn test_high_lightness() {
        let gamut = CieLChGamut::Cie1931D65;
        // Lightness at 99.5, and hue 0.0, should return 2.5 chroma
        // Interpolate between 100.0 and the last value value in the first column of the gamut table
        approx::assert_abs_diff_eq!(gamut.max_chroma(99.5, 0.0).unwrap(), 2.5, epsilon = 1e-6);
        // Lightness at 0.5, and hue 0.5*5 = 2.5, should return 6.5 chroma as well
        // Interpolate between 100 and last value in the first column of the gamut table
        approx::assert_abs_diff_eq!(gamut.max_chroma(99.5, 2.5).unwrap(), 2.5, epsilon = 1e-6);
    }
}
