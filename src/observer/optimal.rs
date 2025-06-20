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

/// `OptimalColors` represents a matrix of tristimulus values derived from optimal spectra,
/// which are theoretical reflectance curves that define the outer boundary (gamut limit)
/// of a color space under a given illuminant and observer.
///
/// The colors are stored in an `NS × NS` matrix, where `NS` is the number of spectral samples.  The
/// values are computed by integrating blocks of the spectral locus, effectively simulating the
/// addition of contiguous spectral bands, with increasing width.  Each row in the matrix
/// corresponds to a different block width, and each column corresponds to a different starting
/// wavelength in the spectrum.  The block wraps around the spectrum to ensure continuity, and to
/// represent the full range of possible colors, including purples.
///
/// The accompanying `XYZ` value stores the reference white used in the computations.
pub struct OptimalColors(XYZ, DMatrix<Vector3<f64>>);

impl Observer {
    /// Computes the optimal colors for the observer under the specified reference white.
    /// The resulting `OptimalColors` contains a matrix of tristimulus values
    /// derived from optimal spectra.
    /// # Arguments
    /// * `ref_white` - The reference white illuminant to use for the calculations.
    /// # Returns
    /// * `OptimalColors` - A struct containing the reference white and the matrix of optimal colors.
    /// # Notes
    /// * The optimal colors are computed by integrating blocks of the spectral locus,
    ///   simulating the addition of contiguous spectral bands with increasing width.
    /// * Each row in the resulting matrix corresponds to a different block width,
    ///   The block wraps around the spectrum to ensure continuity and to represent the full range of possible colors.
    /// * This is a theoretical construct and may not correspond to physically realizable colors.
    /// * As it calculates about 160K colors, it may be computationally intensive.
    pub fn optimal_colors(&self, ref_white: CieIlluminant) -> OptimalColors {
        let mut optcol: DMatrix<Vector3<f64>> =
            DMatrix::from_fn(NS - 1, NS, |_, _| Vector3::zeros());
        let spectral_locus = self.spectral_locus(ref_white);
        let white_point = spectral_locus[0].1.white_point();

        for r in 0..NS - 1 {
            // row
            for c in 0..NS {
                // column, scroll block through spectrum from blue to red
                for w in 0..=r {
                    // integrate block with width r+1
                    let i = (c + w) % NS; // wrap around
                    optcol[(r, c)] += spectral_locus[i].1.xyz().xyz;
                }
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

            if l_u8 > 0 && l_u8 < 100 {
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
}
