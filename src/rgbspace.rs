use std::{collections::HashMap, sync::OnceLock};

use crate::xyz::Chromaticity;
use crate::{
    colorant::Colorant, data::illuminants::D65, data::observers::CIE1931, gamma::GammaCurve,
    illuminant::Illuminant, rgb::gaussian_filtered_primaries, spectrum::Spectrum,
    std_illuminants::StdIlluminant, stimulus::Stimulus,
};
use strum_macros::EnumIter;
use wasm_bindgen::prelude::wasm_bindgen;

// The display P3 red coordinate is outside the CIE 1931 gamut using the CIE 1931 1 nanometer
// dataset as provided by the CIE.  To still match it, it's desatured it adding white. It also mixes
// in a bit of blue (1.4E-6) to get to to the target.
const D: f64 = 0.0005620; // Desaturation ratio
const D65X: f64 = 0.312_738;
const D65Y: f64 = 0.329_052;

#[derive(Debug, Clone, Copy, Default, EnumIter, PartialEq)]
#[wasm_bindgen]
/**
A Light Weight tag, representing an RGB color space.
Used for example in the RGB value set, to identify the color space being used.
 */
pub enum RgbSpace {
    #[default]
    SRGB,
    ADOBE,
    DisplayP3,
}

impl RgbSpace {
    /// Returns the name of the RGB color space.
    pub fn name(&self) -> &'static str {
        match self {
            Self::SRGB => "sRGB",
            Self::ADOBE => "Adobe RGB",
            Self::DisplayP3 => "Display P3",
        }
    }

    /// Obtain reference to the RgbSpace data.
    ///
    /// An `RgbSpace` contains the primary and white spectra, and the gamma
    /// curve function used for encoding and decoding RGB values into a data
    /// stream or image pixel.
    pub fn data(&self) -> &RgbSpaceData {
        match self {
            Self::SRGB => RgbSpaceData::srgb(),
            Self::ADOBE => RgbSpaceData::adobe_rgb(),
            Self::DisplayP3 => RgbSpaceData::display_p3(),
        }
    }

    /// Returns the chromaticity coordinates for the primaries (red, green and blue) of the
    /// RGB colorspace, using the CIE 1931 standard observer.
    pub const fn primaries_chromaticity(&self) -> [Chromaticity; 3] {
        const SRGB_PRIMARIES: [Chromaticity; 3] = [
            Chromaticity::new(0.64, 0.33),
            Chromaticity::new(0.3, 0.6),
            Chromaticity::new(0.15, 0.06),
        ];
        const ADOBE_PRIMARIES: [Chromaticity; 3] = [
            Chromaticity::new(0.64, 0.33),
            Chromaticity::new(0.21, 0.71),
            Chromaticity::new(0.15, 0.06),
        ];
        const DISPLAY_P3_PRIMARIES: [Chromaticity; 3] = [
            Chromaticity::new((1.0 - D) * 0.68 + D * D65X, (1.0 - D) * 0.32 + D * D65Y),
            Chromaticity::new(0.265, 0.69),
            Chromaticity::new(0.15, 0.06),
        ];
        match self {
            Self::SRGB => SRGB_PRIMARIES,
            Self::ADOBE => ADOBE_PRIMARIES,
            Self::DisplayP3 => DISPLAY_P3_PRIMARIES,
        }
    }

    /// Returns the reference white point illuminant for this color space.
    pub fn white(&self) -> StdIlluminant {
        match self {
            Self::SRGB => StdIlluminant::D65,
            Self::ADOBE => StdIlluminant::D65,
            Self::DisplayP3 => StdIlluminant::D65,
        }
    }
}

/**
Spectrally based color space, using spectral representations of the primaries and the
reference white.

Using the CIE 1931 standard observer, using a wavelength domain from 380 top 780
nanometer with 1 nanometer steps, these result in their usual chromaticity
values.  The most common _sRGB_ color space is obtained using the
`RgbSpace::srgb()` constructor. For this instance, the blue and green primaries
are direct Gaussian-filtered D65 spectra. A mixture of the blue primary and a
G1aussian-filtered red component is used for the red primary. Similar
constructors are provided for the `Adobe` and `DisplayP3` color spaces.

The benefit of spectral primaries is that color management and color profiles
can use updated Colorimetric Observers, such as the Cone-Fundamental based CIE
2015 observers, which don't have the CIE 1931 deficiencies. For example, they
can also be optimized for special observers by considering an observer's age or
health conditions.
*/
pub struct RgbSpaceData {
    pub(crate) primaries: [Stimulus; 3],
    pub(crate) white: StdIlluminant,
    pub(crate) gamma: GammaCurve,
}

impl RgbSpaceData {
    /**
    Creates a new `RgbSpace` using an array of three spectra as first
    argument, representing the spectral distributions of the red, green, and
    blue primary stimuli, the white spectrum as second argument, and a gamma
    function as last arguments.  The spectral primaries or white spectrum
    can have arbitrary power, as they are normalized before using in
    calculation, in particular when using the RGB to XYZ transforms, and
    vice versa (see the `rgb2xyz(rgbid: &RgbSpaceId)` and `xyz2rgb(rgbid: &RgbSpaceId)`
    methods of `Observer`).
    */
    pub fn new(primaries: [Stimulus; 3], white: StdIlluminant, gamma: GammaCurve) -> Self {
        Self {
            primaries,
            white,
            gamma,
        }
    }

    /**
     Get primaries as colorants.

     This is the primaries stimulus divided by the reference white.

     # Panics

     If the RGB space reference white has 0.0 values in it, then this cause a division by zero.
    */
    pub fn primaries_as_colorants(&self) -> [Colorant; 3] {
        let white = self
            .white
            .illuminant()
            .clone()
            .set_illuminance(&CIE1931, 100.0)
            .0;
        // RGB primaries defined with reference to CIE1931, and 100 cd/m2.
        let sa = self.primaries.each_ref().map(|v| &v.0 / &white);
        sa.map(Colorant)
    }

    /**
    The sRGB color space, created by HP and Microsoft in 1996.  It is the
    default color space used in an image, or a Web-page, if no color space
    is specified in a tag or in a color profile.
    In this instance, the primaries are composed of Daylight D65 filtered Gaussian filters.
    For the green and blue these are single Gaussians.
    For red primary a combination of a blue and red Gaussian is used, with
    the blue being equal to the blue primary.
    */
    pub fn srgb() -> &'static RgbSpaceData {
        static SRGB: OnceLock<RgbSpaceData> = OnceLock::new();
        /*
        Red Gaussian center wavelength and width, and blue lumen fraction
        added to match the red primary sRGB color space specifiation.
        See the `examples/primaries/main.rs` how these were calculated.
        */
        const RED: [f64; 3] = [627.2101041540204, 23.38636113607498, 0.006839349789397155];

        // Blue and Green Gaussian filters.
        const GREEN: [f64; 2] = [541.2355906535001, 33.66683554608116];
        const BLUE: [f64; 2] = [398.0273721579992, 52.55338039394701];

        SRGB.get_or_init(|| {
            let primaries = gaussian_filtered_primaries(&D65, RED, GREEN, BLUE);
            let white = StdIlluminant::D65;
            let gamma =
                GammaCurve::new(vec![2.4, 1.0 / 1.055, 0.055 / 1.055, 1.0 / 12.92, 0.04045]);
            Self {
                primaries,
                white,
                gamma,
            }
        })
    }
    /**
    The Adobe RGB color space.

    The primaries are composed of Daylight D65 filtered Gaussian filters.
    For the green and blue these are single Gaussians.
    For red primary a combination of a blue and red Gaussian is used, with
    the blue being equal to the blue primary.
    */
    pub fn adobe_rgb() -> &'static RgbSpaceData {
        static ADOBE_RGB: OnceLock<RgbSpaceData> = OnceLock::new();
        // Red Gaussian center wavelength and width, and blue lumen fraction
        // added to match the red primary sRGB color space specifiation.
        // See the `examples/primaries/main.rs` how these were calculated.
        const RED: [f64; 3] = [627.2101041540204, 23.38636113607498, 0.006839349789397155];

        // Blue and Green Gaussian filters.
        const GREEN: [f64; 2] = [531.1505632157933, 20.61013919689458];
        const BLUE: [f64; 2] = [398.0273721579992, 52.55338039394701];

        ADOBE_RGB.get_or_init(|| {
            let primaries = gaussian_filtered_primaries(&D65, RED, GREEN, BLUE);
            let white = StdIlluminant::D65;
            let gamma = GammaCurve::new(vec![563.0 / 256.0]);
            // See https://en.wikipedia.org/wiki/Adobe_RGB_color_space#ICC_PCS_color_image_encoding
            Self {
                primaries,
                white,
                gamma,
            }
        })
    }

    /**
    The Display P3 color space.
    The primaries are composed of Daylight D65 filtered Gaussian filters.
    For the green and blue these are single Gaussians.
    For red primary a combination of a blue and red Gaussian is used, with
    the blue being equal to the blue primary.
    */
    pub fn display_p3() -> &'static RgbSpaceData {
        static DISPLAY_P3: OnceLock<RgbSpaceData> = OnceLock::new();
        // Red Gaussian center wavelength and width, and blue lumen fraction
        // added to match the red primary sRGB color space specifiation.
        // See the `examples/primaries/main.rs` how these were calculated.
        const RED: [f64; 3] = [637.9658554073235, 23.274151215539906, 1.4261447449730065e-6];

        // Blue and Green Gaussian filters.
        const GREEN: [f64; 2] = [539.8416064376327, 21.411199475289777];
        const BLUE: [f64; 2] = [398.0273721579992, 52.55338039394701];

        DISPLAY_P3.get_or_init(|| {
            let primaries = gaussian_filtered_primaries(&D65, RED, GREEN, BLUE);
            let white = StdIlluminant::D65;
            let gamma =
                GammaCurve::new(vec![2.4, 1.0 / 1.055, 0.055 / 1.055, 1.0 / 12.92, 0.04045]);
            Self {
                primaries,
                white,
                gamma,
            }
        })
    }
}

#[cfg(test)]
mod rgbspace_tests {
    use crate::{prelude::*, rgb};
    use approx::assert_ulps_eq;
    use strum::IntoEnumIterator;

    #[test]
    /// Check color points of the primaries, as calculated from the space's
    /// spectra, to the targets returned by the `RgbSpace::primaries_chromaticity()` method.
    fn primaries_chromaticity_match_stimulus_spectrum() {
        for space in RgbSpace::iter() {
            let primaries_chromaticity = space.primaries_chromaticity();
            let primaries_colorants = &space.data().primaries;
            assert_eq!(primaries_chromaticity.len(), primaries_colorants.len());

            let iter = primaries_chromaticity.into_iter().zip(primaries_colorants);
            for (chromaticity, colorant) in iter {
                let computed_chromaticity = CIE1931
                    .xyz_from_spectrum(&colorant.spectrum())
                    .chromaticity();
                assert_ulps_eq!(
                    chromaticity.to_array().as_ref(),
                    computed_chromaticity.to_array().as_ref(),
                    epsilon = 1E-5
                );
            }
        }
    }
}
