//! Provides functionality for handling relative XYZ color space gamuts with different observers and illuminants.
//!
//! This module implements gamut handling for various CIE standard observers and illuminants, allowing for:
//! - Chromaticity binning and conversion
//! - Maximum luminance determination for specific chromaticity coordinates
//! - Conversion between bin coordinates and RelXYZ values
//!
//! # Features
//! - Default support for CIE 1931 2° observer with D50 and D65 illuminants
//! - Optional support for additional observers (CIE 1964, CIE 2015) via "supplemental-observers" feature
//! - High-resolution chromaticity binning (2000x2000 bins)

use std::collections::HashMap;

use crate::illuminant::CieIlluminant;
use crate::observer::Observer;
use crate::traits::Light;
use crate::xyz::RelXYZ;
use std::sync::LazyLock;

/// Configuration constant defining the resolution for chromaticity binning.
/// Higher values provide more precise gamut boundary representation at the cost of memory usage.
/// Resolution for chromaticity binning (2000 bins per axis)
const CHROMATICITY_BINS: u16 = 2000;

/// Maximum valid bin index, derived from CHROMATICITY_BINS.
/// Used for clamping bin values to valid range.
/// Maximum valid bin index (CHROMATICITY_BINS - 1)
const MAX_BIN: u16 = CHROMATICITY_BINS - 1;

/// Scaling factor used for converting between chromaticity values and bin indices.
/// Provides linear mapping between [0,1] chromaticity range and [0,CHROMATICITY_BINS] bin range.
/// Scale factor for converting chromaticity to bin index
const BIN_SCALE: f64 = CHROMATICITY_BINS as f64;

/// Stores gamut data for a specific combination of observer and illuminant.
/// Contains maximum luminance values indexed by chromaticity coordinates.
macro_rules! generate_relxyz_gamuts {
    ($(($identifier:ident, $observer:ident, $illuminant:ident)),* $(,)?) => {
        $(
            static $identifier: LazyLock<RelXYZGamutData> = LazyLock::new(|| {
              //  let opt_colors = Observer::$observer.optimal_colors(CieIlluminant::$illuminant);
              //  let hashmap = opt_colors.max_luminance_per_chromaticity_bin();
                RelXYZGamutData::new(
                    Observer::$observer,
                    CieIlluminant::$illuminant,
                )
            });
        )*
    };
}

generate_relxyz_gamuts!(
    (RELXYZ_GAMUT_CIE1931_D50_DATA, Cie1931, D50),
    (RELXYZ_GAMUT_CIE1931_D65_DATA, Cie1931, D65),
);

#[cfg(feature = "supplemental-observers")]
generate_relxyz_gamuts!(
    (CIEXYZ_GAMUT_CIE1964_D50_DATA, Cie1964, D50),
    (CIEXYZ_GAMUT_CIE1964_D65_DATA, Cie1964, D65),
    (CIEXYZ_GAMUT_CIE2015_D50_DATA, Cie2015, D50),
    (CIEXYZ_GAMUT_CIE2015_D65_DATA, Cie2015, D65),
    (CIEXYZ_GAMUT_CIE2015_10_D50_DATA, Cie2015_10, D50),
    (CIEXYZ_GAMUT_CIE2015_10_D65_DATA, Cie2015_10, D65),
);

struct RelXYZGamutData {
    illuminant: CieIlluminant,
    observer: Observer,
    max_luminances: HashMap<[u16; 2], u16>,
}

impl RelXYZGamutData {
    /*
    const fn new(
        observer: Observer,
        illuminant: CieIlluminant,
        max_luminances: HashMap<[u16; 2], u16>,
    ) -> Self {
        RelXYZGamutData {
            illuminant,
            observer,
            max_luminances,
        }
    }
     */

    /// Computes the maximum luminance (Y) for each chromaticity bin in the CIE xyY color space,
    /// using the optimal color set for this observer and reference white.
    ///
    /// This mapping is useful for:
    /// - Visualizing the outer boundary of the visible gamut in chromaticity diagrams.
    /// - Determining the maximum achievable luminance at each chromaticity point for theoretical (optimal) colors.
    /// - Gamut mapping and boundary visualization tasks.
    ///
    /// # Returns
    /// A `HashMap<[u16; 2], u16>` mapping each chromaticity bin `[x_bin, y_bin]` to the maximum luminance (Y) found in that bin.
    pub fn new(observer: Observer, illuminant: CieIlluminant) -> Self {
        let opt_colors = observer.optimal_colors(illuminant);
        let mut max_luminances = HashMap::new();
        for xyz in opt_colors.colors().iter() {
            let rel_xyz = RelXYZ::from_vec(*xyz, opt_colors.white_point());
            let [xx, yy] = rel_xyz.xyz().chromaticity().to_array();

            // xx < 1.0, and yy<1.0 because they are chromaticity coordinates.
            let x_bin = Self::chromaticity_to_bin(xx);
            let y_bin = Self::chromaticity_to_bin(yy);

            // Scale the luminance to a u16 value
            // The luminance is scaled from [0, 100.0] to [0, 65535]
            let y_max_u16 = Self::luminance_to_l_bin(rel_xyz.xyz().y());

            // Only add if the discrete chromaticity coordinates and luminance are valid.
            // This does a to and from CieLab round trip, which is not ideal,
            // but it is a reasonable approximation for the purpose of this mapping.
            let rxyz_for_bin =
                Self::bins_to_rel_xyz_static(observer, illuminant, x_bin, y_bin, y_max_u16);
            if rxyz_for_bin.is_valid() {
                // Insert or update the maximum luminance for this chromaticity bin
                max_luminances
                    .entry([x_bin, y_bin])
                    .and_modify(|existing| {
                        if *existing < y_max_u16 {
                            *existing = y_max_u16;
                        }
                    })
                    .or_insert(y_max_u16);
            }
        }
        RelXYZGamutData {
            illuminant,
            observer,
            max_luminances,
        }
    }

    pub fn bins_to_rel_xyz_static(
        observer: Observer,
        illuminant: CieIlluminant,
        x_bin: u16,
        y_bin: u16,
        y_max_u16: u16,
    ) -> RelXYZ {
        let x = Self::bin_to_chromaticity(x_bin);
        let y = Self::bin_to_chromaticity(y_bin);
        let z = 1.0 - x - y;
        let l = y_max_u16 as f64 / u16::MAX as f64 * 100.0;
        let scale = l / y.max(f64::EPSILON);
        RelXYZ::new([x * scale, l, z * scale], illuminant.xyzn100(observer))
    }

    pub fn chromaticity_to_bin(value: f64) -> u16 {
        (value.clamp(0.0, 1.0) * BIN_SCALE).floor() as u16
    }

    pub fn bin_to_chromaticity(bin: u16) -> f64 {
        (bin.min(MAX_BIN) as f64) / BIN_SCALE
    }

    pub fn luminance_to_l_bin(luminance: f64) -> u16 {
        // Convert luminance to a u16 value
        (luminance / 100.0 * u16::MAX as f64).round() as u16
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, strum::Display)]
/// Represents available RelXYZ gamut configurations for different observer/illuminant combinations.
///
/// Provides methods for:
/// - Accessing observer and illuminant information
/// - Converting between chromaticity values and bin indices
/// - Retrieving maximum luminance values for specific chromaticity coordinates
/// - Converting bin coordinates to RelXYZ colors
///
/// The default configuration is CIE 1931 2° observer with D65 illuminant.
/// Additional configurations are available with the "supplemental-observers" feature.
pub enum RelXYZGamut {
    #[default]
    Cie1931D65,
    Cie1931D50,
    #[cfg(feature = "supplemental-observers")]
    Cie1964D50,
    #[cfg(feature = "supplemental-observers")]
    Cie1964D65,
    #[cfg(feature = "supplemental-observers")]
    Cie2015D50,
    #[cfg(feature = "supplemental-observers")]
    Cie2015D65,
    #[cfg(feature = "supplemental-observers")]
    Cie2015_10D50,
    #[cfg(feature = "supplemental-observers")]
    Cie2015_10D65,
}

impl RelXYZGamut {
    pub fn observer(&self) -> Observer {
        self.data().observer
    }

    pub fn illuminant(&self) -> CieIlluminant {
        self.data().illuminant
    }

    pub fn bins_to_rel_xyz(&self, x_bin: u16, y_bin: u16, y_max_u16: u16) -> RelXYZ {
        RelXYZGamutData::bins_to_rel_xyz_static(
            self.observer(),
            self.illuminant(),
            x_bin,
            y_bin,
            y_max_u16,
        )
    }

    fn data(&self) -> &RelXYZGamutData {
        match self {
            RelXYZGamut::Cie1931D50 => &RELXYZ_GAMUT_CIE1931_D50_DATA,
            RelXYZGamut::Cie1931D65 => &RELXYZ_GAMUT_CIE1931_D65_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie1964D50 => &CIEXYZ_GAMUT_CIE1964_D50_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie1964D65 => &CIEXYZ_GAMUT_CIE1964_D65_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie2015D50 => &CIEXYZ_GAMUT_CIE2015_D50_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie2015D65 => &CIEXYZ_GAMUT_CIE2015_D65_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie2015_10D50 => &CIEXYZ_GAMUT_CIE2015_10_D50_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie2015_10D65 => &CIEXYZ_GAMUT_CIE2015_10_D65_DATA,
        }
    }

    pub fn max_luminance(&self, x: u16, y: u16) -> Option<u16> {
        self.data().max_luminances.get(&[x, y]).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_relxyz_gamut_data() {
        let gamut = RelXYZGamut::Cie1931D65;
        assert_eq!(gamut.observer(), Observer::Cie1931);
        assert!(!gamut.data().max_luminances.is_empty());
    }

    #[test]
    fn test_max_luminance_for_bin() {
        let gamut = RelXYZGamut::Cie1931D65;
        for x in 200..=300 {
            for y in 200..=300 {
                let max_luminance = gamut.max_luminance(x, y);
                if let Some(luminance) = max_luminance {
                    println!("Max luminance found for bin ({}, {}): {}", x, y, luminance);
                    return;
                }
            }
        }
        panic!("No max luminance found for any bin in the range 200-300");
    }

    #[test]
    fn test_chromaticity_to_bin() {
        assert_eq!(RelXYZGamutData::chromaticity_to_bin(0.0), 0);
        assert_eq!(RelXYZGamutData::chromaticity_to_bin(0.5), 1000);
        assert_eq!(RelXYZGamutData::chromaticity_to_bin(0.9995), 1999);

        // Test edge cases
        assert_eq!(RelXYZGamutData::chromaticity_to_bin(-0.1), 0);
        assert_eq!(RelXYZGamutData::chromaticity_to_bin(1.1), 2000);

        // Test roundtrip
        let test_values = [0.1, 0.25, 0.5, 0.75, 0.9];
        for &x in &test_values {
            let bin = RelXYZGamutData::chromaticity_to_bin(x);
            let x_restored = RelXYZGamutData::bin_to_chromaticity(bin);
            assert_abs_diff_eq!(x, x_restored, epsilon = 0.0005);
        }
    }

    #[test]
    fn test_bins_to_rel_xyz() {
        let white_point = Observer::Cie1931.xyz_d65();
        let gamut = RelXYZGamut::Cie1931D65;

        // Test center point (0.3, 0.5)
        let [xx, yy] = white_point.chromaticity().to_array();
        let xbin = RelXYZGamutData::chromaticity_to_bin(xx);
        let ybin = RelXYZGamutData::chromaticity_to_bin(yy);
        let center = gamut.bins_to_rel_xyz(xbin, ybin, 32768);
        assert_abs_diff_eq!(center.xyz().x(), 47.49, epsilon = 0.005);
        assert_abs_diff_eq!(center.xyz().y(), 50.0, epsilon = 0.005);
        assert_abs_diff_eq!(center.xyz().z(), 54.48, epsilon = 0.005);
    }

    #[test]
    fn test_max_luminance() {
        let max_luminances = RelXYZGamut::Cie1931D65;

        let [x, y] = max_luminances
            .observer()
            .xyz_d65()
            .chromaticity()
            .to_array();
        let x_bin = (x * 2000.0).round() as u16;
        let y_bin = (y * 2000.0).round() as u16;
        let y_max_u16 = max_luminances.max_luminance(x_bin, y_bin).unwrap();
        assert!(
            y_max_u16 == u16::MAX,
            "Expected maximum luminance of 255 for D65 reference white"
        );
    }
}
