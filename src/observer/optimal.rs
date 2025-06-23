//! Optimal color computations for colorimetry observers.
//!
//! This module provides functionality to compute the theoretical gamut boundary
//! (optimal colors) for a given observer and illuminant. It defines the
//! [`OptimalColors`] struct, which stores a matrix of XYZ tristimulus values
//! derived from optimal spectra, as well as methods for generating and working
//! with these values.
//!
//! The optimal color computations are based on the observer’s color matching
//! functions and a specified reference white, and are useful for visualizing
//! the limits of color spaces under different viewing conditions.

use nalgebra::{DMatrix, Vector3};
use std::collections::HashMap;

use super::Observer;
use crate::{
    lab::CieLab,
    prelude::CieIlluminant,
    spectrum::NS,
    xyz::{RelXYZ, XYZ},
};

/// Represents the set of optimal colors in colorimetry.
///
/// Optimal colors define the theoretical limits of color saturation (chroma) for given
/// lightness and hue values. They are based on binary reflectance spectra that reflect
/// either 0% or 100% of light at each wavelength, representing the most vivid colors
/// physically possible under the constraints of human vision.
///
/// These colors form the boundary of the CIE color spaces (such as CIE 1931 or CIELAB),
/// and are used to compute:
///
/// - The spectral locus in chromaticity diagrams.
/// - The maximum chroma values at given lightness and hue levels in perceptual spaces like CIELAB.
///
/// Although optimal colors do not correspond to real-world pigments or displays, they are
/// essential for defining the theoretical gamut of visible colors.
pub struct OptimalColors(XYZ, DMatrix<Vector3<f64>>);

impl Observer {
    /// Computes the optimal colors for the observer under the specified reference white.
    ///
    /// The resulting `OptimalColors` contains a matrix of tristimulus values
    /// derived from optimal spectra.
    ///
    /// # Arguments
    /// * `ref_white` - The reference white illuminant to use for the calculations.
    ///
    /// # Returns
    /// * `OptimalColors` - A struct containing the reference white and the matrix of optimal colors.
    ///
    /// # Notes
    /// * The optimal colors are computed by applying a rectangular filter to the spectral locus,
    ///   simulating the addition of contiguous spectral bands with increasing width.
    /// * Each row in the resulting matrix corresponds to a different filter width,
    ///   and the filter wraps around the spectrum to ensure continuity and to represent the full range of possible colors.
    /// * The number of spectral samples, `NS`, is fixed at 401 throughout the library.
    ///   - The resulting matrix has shape `(NS-1, NS)` = (400, 401).
    ///   - Each column corresponds to a starting wavelength, and each row to a filter width.
    /// * This is a theoretical construct and may not correspond to physically realizable colors.
    /// * As it calculates about 160,000 colors, it may be computationally intensive.
    pub fn optimal_colors(&self, ref_white: CieIlluminant) -> OptimalColors {
        let mut optcol: DMatrix<Vector3<f64>> = DMatrix::zeros(NS - 1, NS);
        let spectral_locus = self.spectral_locus(ref_white);
        let white_point = spectral_locus[0].1.white_point();

        for c in 0..NS {
            optcol[(0, c)] = spectral_locus[c].1.xyz().xyz; // first row is the spectral locus
        }

        for r in 1..NS - 1 {
            for c in 0..NS {
                let i = (c + r) % NS; // wrap around
                optcol[(r, c)] = optcol[(r - 1, c)] + optcol[(0, i)];
            }
        }
        OptimalColors(white_point, optcol)
    }
}

impl OptimalColors {
    pub const LAB_MAX_CHROME_LEN: usize = 7128;
    pub const NH: usize = 72;
    pub const NL: usize = 99;

    pub fn white_point(&self) -> XYZ {
        self.0
    }

    pub fn colors(&self) -> &DMatrix<Vector3<f64>> {
        &self.1
    }

    pub fn observer(&self) -> Observer {
        self.0.observer
    }

    /// Computes the maximum chroma values across lightness–hue combinations in the CIE LCh color space.
    ///
    /// Returns a `HashMap` where each key is a `(L, H)` pair:
    /// - `L` is lightness (rounded to the nearest integer in the range 0–100),
    /// - `H` is the hue category, with the hue angles in 72 bins, each covering 5 degrees.
    ///
    /// Each value in the map is the highest chroma (`C`) observed for that `(L, H)` bin,
    /// based on the XYZ values stored in the `OptimalColors` matrix, and dimmer values to fill the bins.
    ///
    /// This method is useful for estimating the outer boundary (gamut limit) of the LCh color
    /// space, especially when visualizing or constraining perceptually uniform color
    /// representations.
    pub fn cielab_max_chromas(&self) -> HashMap<(u8, u8), u8> {
        let mut map = HashMap::new();
        let white_point = self.0;

        let mut add_xyz = |xyz: Vector3<f64>, scale: f64| {
            let xyz_scaled = xyz * scale;
            let rel_xyz = RelXYZ::from_vec(xyz_scaled, white_point);
            let lab = CieLab::from_xyz(rel_xyz);
            let [l, c, h] = lab.lch();

            let l_u8 = l.round() as u8;
            let h_u8 = (h / 5.0).floor() as u8;
            let c_u8 = c.floor() as u8;

            // Only add to the map if lightness is in the range 1–99
            // and hue is in the range 0–71 (72 bins covering 360 degrees).
            // Lightness values 0 and 100 are excluded, as they represent edge cases
            // that are always considered in gamut for all hues, and with zero chroma.
            if (1..100).contains(&l_u8) && h_u8 < Self::NH as u8 {
                map.entry((l_u8, h_u8))
                    .and_modify(|existing| {
                        if *existing < c_u8 {
                            *existing = c_u8;
                        }
                    })
                    .or_insert(c_u8);
            }
        };

        // fill the map with the original colors, and then with dimmer versions to cover the full gamut
        for xyz in self.1.iter() {
            let mut scale = 1.0;
            for _ in 1..=7 {
                add_xyz(*xyz, scale);
                scale /= 1.61803; // golden ratio
            }
        }
        map
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observer::Observer::Cie1931;
    use approx::assert_abs_diff_eq;
    use strum::IntoEnumIterator;

    #[test]
    fn test_optimum_colors() {
        let observer = Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = observer.optimal_colors(ref_white);

        let first: Vec<XYZ> = observer
            .spectral_locus(ref_white)
            .iter()
            .map(|(_w, xyz)| xyz.xyz())
            .collect();

        first
            .iter()
            .map(|xyz| xyz.xyz)
            .zip(opt_colors.colors().row(0).iter())
            .for_each(|(x, y)| {
                assert_abs_diff_eq!(x, y, epsilon = 1e-10);
            });

        dbg!(
            opt_colors.colors().row(NS - 2)[0],
            opt_colors.colors().row(NS - 2)[NS - 1]
        );
    }

    #[test]
    fn test_cielch() {
        let observer = Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = observer.optimal_colors(ref_white);
        let cielch = opt_colors.cielab_max_chromas();
        for l in 1..=99 {
            for h in 0..72 {
                if !cielch.contains_key(&(l, h)) {
                    eprintln!("bin {}, {} not populated", l, h);
                }
            }
        }
        assert_eq!(cielch.len(), OptimalColors::LAB_MAX_CHROME_LEN);
    }

    #[test]
    #[ignore = "This test is computationally intensive and may take a long time to run."]
    fn test_cielch_all_observers() {
        for ref_white in [CieIlluminant::D65, CieIlluminant::D50] {
            for obs in Observer::iter() {
                let opt_colors = obs.optimal_colors(ref_white);
                let cielch = opt_colors.cielab_max_chromas();
                assert_eq!(
                    cielch.len(),
                    OptimalColors::LAB_MAX_CHROME_LEN,
                    "Observer: {}, Reference White: {:?} has incorrect cielch length",
                    obs.as_ref(),
                    ref_white
                );
            }
        }
    }

    #[test]
    fn test_values() {
        let observer = Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = observer.optimal_colors(ref_white);
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
        let cielab_max_chromas = opt_colors.cielab_max_chromas();
        for &[l, h, expected] in WANTS {
            let &c = cielab_max_chromas.get(&(l, h)).unwrap_or(&0);
            assert_eq!(
                c, expected,
                "Expected chroma for (L={}, H={}) to be {}, but got {}",
                l, h, expected, c
            );
        }
    }

    #[test]
    fn test_optimal_colors_matrix_shape_and_first_row() {
        let observer = Observer::Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = observer.optimal_colors(ref_white);
        let matrix = opt_colors.colors();

        // Check matrix shape
        assert_eq!(matrix.nrows(), NS - 1, "Matrix should have NS-1 rows");
        assert_eq!(matrix.ncols(), NS, "Matrix should have NS columns");

        // Check first row matches spectral locus
        let spectral_locus = observer.spectral_locus(ref_white);
        for (c, (_, sl_val)) in spectral_locus.iter().enumerate() {
            let expected = sl_val.xyz().xyz;
            let actual = matrix.row(0)[c];
            assert_abs_diff_eq!(expected, actual, epsilon = 1e-10);
        }
    }
}
