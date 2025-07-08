use std::fmt;

use crate::xyz::{Chromaticity, XYZ};
use crate::{
    colorant::Colorant,
    illuminant::CieIlluminant,
    math::Triangle,
    observer::Observer::{self, Cie1931},
    rgb::gamma::GammaCurve,
    stimulus::Stimulus,
};
use strum::{AsRefStr, EnumCount, EnumIter};

mod srgb;
use srgb::SRGB;

mod display_p3;
use display_p3::DISPLAY_P3;

mod adobe_rgb;
use adobe_rgb::ADOBE_RGB;

mod cie_rgb;
use cie_rgb::CIE_RGB;

/// Spectrally based color space, using spectral representations of the primaries and the
/// reference white.
///
/// Using the CIE 1931 standard observer, using a wavelength domain from 380 top 780
/// nanometer with 1 nanometer steps, these result in their usual chromaticity
/// values.  The most common _sRGB_ color space is obtained using the
/// `RgbSpace::srgb()` constructor. For this instance, the blue and green primaries
/// are direct Gaussian-filtered D65 spectra. A mixture of the blue primary and a
/// G1aussian-filtered red component is used for the red primary. Similar
/// constructors are provided for the `Adobe` and `DisplayP3` color spaces.
///
/// The benefit of spectral primaries is that color management and color profiles
/// can use updated Colorimetric Observers, such as the Cone-Fundamental based CIE
/// 2015 observers, which don't have the CIE 1931 deficiencies. For example, they
/// can also be optimized for special observers by considering an observer's age or
/// health conditions.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Debug, Clone, Copy, Default, Eq, PartialEq, Hash, EnumIter, AsRefStr)]
#[non_exhaustive]
pub enum RgbSpace {
    #[default]
    SRGB,
    Adobe,
    DisplayP3,
    CieRGB,
}

impl RgbSpace {
    /// Obtain reference to the RgbSpace data.
    ///
    /// An `RgbSpace` contains the primary and white spectra, and the gamma
    /// curve function used for encoding and decoding RGB values into a data
    /// stream or image pixel.
    fn data(&self) -> &RgbSpaceData {
        match self {
            Self::SRGB => &SRGB,
            Self::Adobe => &ADOBE_RGB,
            Self::DisplayP3 => &DISPLAY_P3,
            Self::CieRGB => &CIE_RGB,
        }
    }

    pub fn primaries(&self) -> &[Stimulus; 3] {
        &self.data().primaries
    }

    pub fn name(&self) -> &'static str {
        self.data().name
    }

    /// Returns the reference white point illuminant for this color space.
    pub fn white(&self) -> CieIlluminant {
        self.data().white
    }

    pub fn gamma(&self) -> &GammaCurve {
        &self.data().gamma
    }

    pub fn contains(&self, observer: Observer, chromaticity: &Chromaticity) -> bool {
        let chromaticities = self.chromaticities(observer);
        let triangle = Triangle::new(
            chromaticities[0].to_array(),
            chromaticities[1].to_array(),
            chromaticities[2].to_array(),
        )
        .unwrap(); // validated upon construction
        triangle.contains(chromaticity.x(), chromaticity.y())
    }

    pub fn chromaticities(&self, observer: Observer) -> [Chromaticity; 3] {
        let [[rx, ry], [gx, gy], [bx, by]] = self.data().chromaticities[observer as usize];
        [
            Chromaticity::new(rx, ry),
            Chromaticity::new(gx, gy),
            Chromaticity::new(bx, by),
        ]
    }

    pub fn target_chromaticities_cie1931(&self) -> [Chromaticity; 3] {
        let [[rx, ry], [gx, gy], [bx, by]] = self.data().target_chromaticities_cie1931;
        [
            Chromaticity::new(rx, ry),
            Chromaticity::new(gx, gy),
            Chromaticity::new(bx, by),
        ]
    }

    pub fn white_point(&self, observer: Observer) -> XYZ {
        XYZ::new(self.data().white_points[observer as usize], observer)
    }

    /// Get primaries as colorants.
    ///
    ///  This is the primaries stimulus divided by the reference white.
    ///
    ///  # Panics
    ///
    ///  If the RGB space reference white has 0.0 values in it, then this cause a division by zero.
    pub fn primaries_as_colorants(&self) -> [Colorant; 3] {
        let white = self
            .white()
            .illuminant()
            .clone()
            .set_illuminance(Cie1931, 100.0)
            .0;
        // RGB primaries defined with reference to CIE1931, and 100 cd/m2.
        let sa = self.primaries().each_ref().map(|v| &v.0 / &white);
        sa.map(Colorant)
    }
}

// Implementing the Display trait for RgbSpace to provide a string representation.
// To get the EnumIter to work with strum, use the AsRefStr trait, which is Derived
impl fmt::Display for RgbSpace {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.data().name)
    }
}

struct RgbSpaceData {
    primaries: [Stimulus; 3],
    white: CieIlluminant,
    gamma: GammaCurve,
    name: &'static str,
    target_chromaticities_cie1931: [[f64; 2]; 3],
    chromaticities: [[[f64; 2]; 3]; Observer::COUNT], // Use Observer::COUNT from strum::EnumCount
    white_points: [[f64; 3]; Observer::COUNT],
}

#[cfg(test)]
mod rgbspace_tests {
    use crate::prelude::*;
    use approx::assert_ulps_eq;
    use strum::IntoEnumIterator;
    use Observer::Cie1931;

    #[test]
    /// Check color points of the primaries, as calculated from the space's
    /// spectra, to the targets returned by the `RgbSpace::primaries_chromaticity()` method.
    fn primaries_chromaticity_match_stimulus_spectrum() {
        for space in super::RgbSpace::iter() {
            let primaries_chromaticity = space.target_chromaticities_cie1931();
            let spectral_primaries = &space.data().primaries;
            assert_eq!(primaries_chromaticity.len(), spectral_primaries.len());

            let iter = primaries_chromaticity.into_iter().zip(spectral_primaries);
            for (chromaticity, primary) in iter {
                let computed_chromaticity = Cie1931
                    .xyz_from_spectrum(&primary.spectrum())
                    .chromaticity();
                assert_ulps_eq!(
                    chromaticity.to_array().as_ref(),
                    computed_chromaticity.to_array().as_ref(),
                    epsilon = 5E-4
                );
            }
        }
    }

    #[test]
    fn test_rgb_as_ref_str() {
        assert_eq!(RgbSpace::SRGB.as_ref(), "SRGB");
        assert_eq!(RgbSpace::DisplayP3.as_ref(), "DisplayP3");
        assert_eq!(RgbSpace::Adobe.as_ref(), "Adobe");
        assert_eq!(RgbSpace::CieRGB.as_ref(), "CieRGB");
    }
}
