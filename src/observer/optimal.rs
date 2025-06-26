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

use super::Observer;
use crate::{illuminant::CieIlluminant, spectrum::NS, xyz::XYZ};

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

impl OptimalColors {
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
    pub fn new(observer: Observer, ref_white: CieIlluminant) -> OptimalColors {
        let mut optcol: DMatrix<Vector3<f64>> = DMatrix::zeros(NS - 1, NS);
        let spectral_locus = observer.spectral_locus(ref_white);
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

    pub fn white_point(&self) -> XYZ {
        self.0
    }

    pub fn colors(&self) -> &DMatrix<Vector3<f64>> {
        &self.1
    }

    pub fn observer(&self) -> Observer {
        self.0.observer
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::illuminant::CieIlluminant;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_optimal_colors_matrix_shape_and_first_row() {
        let observer = Observer::Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = OptimalColors::new(observer, ref_white);
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
