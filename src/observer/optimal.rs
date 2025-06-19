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

use std::collections::HashMap;
use nalgebra::{DMatrix, Vector3};

use super::Observer;
use crate::{lab::CieLab, prelude::CieIlluminant, spectrum::NS, xyz::{RelXYZ, XYZ}};

/// `OptimalColors` represents a matrix of tristimulus values derived from optimal spectra,
/// which are theoretical reflectance curves that define the outer boundary (gamut limit)
/// of a color space under a given illuminant and observer.
///
/// The colors are stored in an `NS × NS` matrix, where `NS` is the number of spectral samples.
/// Each entry is a `Vector3<f64>` containing XYZ values. The upper-left triangle contains
/// forward optimal colors (monotonically increasing or decreasing reflectance), while
/// the lower-right triangle contains their inverse counterparts.
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
    pub fn optimal_colors(&self, ref_white: CieIlluminant) -> OptimalColors {
        let mut optcol: DMatrix<Vector3<f64>> = DMatrix::from_fn(NS, NS, |_, _| Vector3::zeros());
        let spectral_locus = self.spectral_locus(ref_white);
        let white_point = spectral_locus[0].1.white_point();
        let white_point_vec = white_point.xyz;

        // fill top left, including diagonal, with direct colors
        for r in 0..NS {
            for c in 0..NS - r {
                for w in 0..=r {
                    optcol[(r, c)] += spectral_locus[c + w].1.xyz().xyz;
                }
            }
        }

        // fill bottom right triangle (inverse colors)
        for r in 0..NS {
            for c in 0..NS {
                if r + c <= NS - 1 {
                    let rr = NS - 1 - r;
                    let cc = NS - 1 - c;
                    optcol[(rr, cc)] = white_point_vec - optcol[(r, c)];
                }
            }
        }
        OptimalColors(white_point, optcol)
    }
}

impl OptimalColors {
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
    /// - `H` is the hue angle (divided by 2 and rounded, giving values in 0–180).
    ///
    /// Each value in the map is the highest chroma (`C`) observed for that `(L, H)` bin,
    /// based on the XYZ values stored in the `OptimalColors` matrix.
    ///
    /// This method is useful for estimating the outer boundary (gamut limit)
    /// of the LCh color space, especially when visualizing or constraining
    /// perceptually uniform color representations.
    pub fn cielab_max_chromas(&self) -> HashMap<(u8, u8), f32> {
        let mut map = HashMap::new();
        let white_point = self.0;

        for xyz in self.1.iter() {
            let rel_xyz = RelXYZ::from_vec(*xyz, white_point);
            let lab = CieLab::from_xyz(rel_xyz);
            let [l, chroma, h] = lab.lch();

            let l = l.round() as u8;
            let h = (h / 2.0).round() as u8;
            let c_f32 = chroma as f32;

            map.entry((l, h))
                .and_modify(|existing| {
                    if *existing < c_f32 {
                        *existing = c_f32;
                    }
                })
                .or_insert(c_f32);
        }
        map
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::observer::Observer::Cie1931;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_optimum_colors() {
        let observer = Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = observer.optimal_colors(ref_white);
        let xyz_ref_white = observer.xyz_d65();

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
        
        let last: Vec<XYZ> = observer
            .spectral_locus(ref_white)
            .iter()
            .rev()
            .map(|(_w, xyz)| xyz_ref_white - xyz.xyz())
            .collect();

        last 
            .iter()
            .map(|xyz| xyz.xyz)
            .zip(opt_colors.colors().row(NS-1).iter())
            .for_each(|(x, y)| {
                println!("{:.5} {:.5} {:.5} {:.5} {:.5} {:.5}", x.x, x.y, x.z, y.x, y.y, y.z);
                assert_abs_diff_eq!(x, y, epsilon = 1e-10);
            });

    }

    #[test]
    fn test_cielch() {
        let observer = Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = observer.optimal_colors(ref_white);
        let cielch = opt_colors.cielab_max_chromas();
        for l in 0..=100 {
            for h in 0..=180 {
                if let Some(c) = cielch.get(&(l, h)) {
                    println!("{}, {}, {}", l, h, c);
                }
            }
        }   
    }
}



