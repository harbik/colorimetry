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
use crate::{
    illuminant::CieIlluminant,
    spectrum::NS,
    xyz::XYZ,
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

   // const RESOLUTION: f64 = 0.0005; // resolution for chromaticity coordinates

    pub fn white_point(&self) -> XYZ {
        self.0
    }

    pub fn colors(&self) -> &DMatrix<Vector3<f64>> {
        &self.1
    }

    pub fn observer(&self) -> Observer {
        self.0.observer
    }

    /*
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

        let max_lums = self.max_luminance_per_chromaticity_bin();
        for (&[x_bin, y_bin], &y_max_u16) in max_lums.iter() {
            let rel_xyz = self.bins_to_rel_xyz(x_bin, y_bin, y_max_u16);
            let [l, c, h] = CieLab::from_xyz(rel_xyz).lch();

            // Round lightness to the nearest integer in the range 0–100
            let l_bin = l.round() as u8;
            // Hue is in degrees, convert to 72 bins (5 degrees each)
            let h_bin = (h / 5.0).round() as u8;
            let c_bin = c.round() as u8;
            // Ensure lightness is within the valid range
            if l_bin <= 100 && h_bin < Self::NH as u8 {
                // Update the maximum chroma for this (L, H) pair
                map.entry((l_bin, h_bin))
                    .and_modify(|c: &mut u8| {
                        if c_bin > *c {
                            *c = c_bin;
                        }
                    })
                    .or_insert(c_bin);
            }
        }
        map
    }

    fn chromaticity_to_xy_bin(value: f64) -> u16 {
        // Scale the chromaticity coordinate to a bin index
        (value.clamp(0.0, 1.0) / Self::RESOLUTION).round() as u16
    }

    fn xy_bin_to_chromaticity(bin: u16) -> f64 {
        // Convert a bin index back to a chromaticity coordinate
        bin as f64 * Self::RESOLUTION
    }

    fn l_xy_bin_to_luminance(y_max_u16: u16) -> f64 {
        y_max_u16 as f64 / u16::MAX as f64 * 100.0
    }

    fn luminance_to_l_xy_bin(luminance: f64) -> u16 {
        // Convert luminance to a u16 value
        (luminance / 100.0 * u16::MAX as f64).round() as u16
    }

    fn bins_to_rel_xyz(&self, x_bin: u16, y_bin: u16, y_max_u16: u16) -> RelXYZ {
        // Convert bins to relative XYZ coordinates
        let x = Self::xy_bin_to_chromaticity(x_bin);
        let y = Self::xy_bin_to_chromaticity(y_bin);
        let z: f64 = 1.0 - x - y; // CIE XYZ constraint
        let l = Self::l_xy_bin_to_luminance(y_max_u16);
        let scale = l / y;
        RelXYZ::new([x * scale, y * scale, z * scale], self.white_point())
    }

    fn lch_bin_to_chromaticity(&self, l_bin: u8, c_bin: u8, h_bin: u8) -> Option<[f64; 2]> {
        // Convert lightness and hue bins to chromaticity coordinates
        let l = l_bin as f64;
        let h = h_bin as f64 * 5.0; // 5 degrees per bin
        let c = c_bin as f64;
        let lab = CieLab::from_lch([l, c, h], self.white_point());
        let xyz = lab.xyz();
        if xyz.is_valid() {
            Some(xyz.xyz().chromaticity().to_array())
        } else {
            // If the XYZ is not valid, return a default value
            None
        }
    }
     */
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::illuminant::CieIlluminant;
    use approx::assert_abs_diff_eq;



    /*
    #[test]
    #[ignore = "This test has output which can be checked."]
    fn test_cielch_hashmap_chromaticity() {
        let observer = Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = observer.optimal_colors(ref_white);
        let cielch = opt_colors.cielab_max_chromas();
        for h in 0..72 {
            print!("[");
            for l in 1..=100 {
                if cielch.contains_key(&(l, h)) {
                    let c = cielch[&(l, h)];
                    if let Some([x, y]) = opt_colors.lch_bin_to_chromaticity(l, c, h) {
                        print!("[{:.5}, {:.5}],", x, y);
                    };
                }
            }
            println!("],");
        }
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
     */

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

    /*
    #[test]
    fn test_chromaticity_to_bin() {
        // Test exact values
        assert_eq!(OptimalColors::chromaticity_to_xy_bin(0.0), 0);
        assert_eq!(OptimalColors::chromaticity_to_xy_bin(0.5), 1000);
        assert_eq!(OptimalColors::chromaticity_to_xy_bin(1.0), 2000);

        // Test edge cases
        assert_eq!(OptimalColors::chromaticity_to_xy_bin(0.0001), 0);
        assert_eq!(OptimalColors::chromaticity_to_xy_bin(0.999999), 2000);

        // Test rounding behavior
        assert_eq!(OptimalColors::chromaticity_to_xy_bin(0.4995), 999);
        assert_eq!(OptimalColors::chromaticity_to_xy_bin(0.5005), 1001);

        // Test roundtrip (bin -> chromaticity -> bin)
        let test_values = [0.1, 0.25, 0.5, 0.75, 0.9];
        for &x in &test_values {
            let bin = OptimalColors::chromaticity_to_xy_bin(x);
            let x_restored = OptimalColors::xy_bin_to_chromaticity(bin);
            assert_abs_diff_eq!(x, x_restored, epsilon = 0.0005);
        }

        // Test invalid values (should clamp)
        assert_eq!(OptimalColors::chromaticity_to_xy_bin(-0.1), 0);
        assert_eq!(OptimalColors::chromaticity_to_xy_bin(1.1), 2000);
    }

    #[test]
    fn test_bins_to_rel_xyz() {
        // Setup test data
        let opt_colors = Cie1931.optimal_colors(CieIlluminant::D65);

        let whitepoint = opt_colors.white_point();
        assert_abs_diff_eq!(whitepoint.x(), 95.042, epsilon = 0.001);
        assert_abs_diff_eq!(whitepoint.y(), 100.000, epsilon = 0.001);
        assert_abs_diff_eq!(whitepoint.z(), 108.861, epsilon = 0.001);

        // Test equal energy (0.3333, 0.3333), with Y = 50.0.
        // bin is 0.3333 * 2000 = 667
        let center = opt_colors.bins_to_rel_xyz(667, 667, 32768);
        assert_abs_diff_eq!(center.xyz().x(), 50.0, epsilon = 0.1); // discretization inaccuracy
        assert_abs_diff_eq!(center.xyz().y(), 50.0, epsilon = 0.1);
        assert_abs_diff_eq!(center.xyz().z(), 50.0, epsilon = 0.1);

        // Test round trip
        let test_bins = [(250, 250), (500, 500), (750, 750)];
        for &(x_bin, y_bin) in &test_bins {
            let rel_xyz = opt_colors.bins_to_rel_xyz(x_bin, y_bin, 32768);
            let [xx, yy] = rel_xyz.xyz().chromaticity().to_array();
            let x_bin_back = OptimalColors::chromaticity_to_xy_bin(xx);
            let y_bin_back = OptimalColors::chromaticity_to_xy_bin(yy);
            assert_eq!(x_bin, x_bin_back);
            assert_eq!(y_bin, y_bin_back);
        }

        // Test luminance scaling
        let d65 = whitepoint.chromaticity();
        let x_bin = OptimalColors::chromaticity_to_xy_bin(d65.x());
        let y_bin = OptimalColors::chromaticity_to_xy_bin(d65.y());
        dbg!(x_bin, y_bin, d65);
        let test_lum = opt_colors.bins_to_rel_xyz(x_bin, y_bin, u16::MAX);
        assert_abs_diff_eq!(whitepoint, test_lum.xyz(), epsilon = 0.2); // discretization inaccuracy
    }
*/
}
