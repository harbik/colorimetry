use crate::{
    colorant::Colorant,
    data::observers::CIE1931,
    illuminant::Illuminant,
    observer::Observer,
    rgbspace::RgbSpace,
    spectrum::Spectrum,
    stimulus::Stimulus,
    traits::{Filter, Light},
    xyz::XYZ,
};
use approx::AbsDiffEq;
use nalgebra::{Matrix3, Vector3};
use std::{borrow::Cow, sync::OnceLock};
use wasm_bindgen::prelude::wasm_bindgen;

/// Representation of a color stimulus in a set of Red, Green, and Blue (RGB) values,
/// representing its relative composition using standard primaries.
///
/// RGB values are commonly used in digital images, with the relative intensity
/// of the primaries defined as three 8-bit values, with range from 0 to 255.
/// As ooposed to CIE XYZ tristimulus values, which used imaginary primaries,
/// displays use real primaries, typically defined in the CIE 1931 diagram.
/// They cover a triangular area, referred to the _color gamut_ of a display.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RGB {
    /// The RGB color space the color values are using. Often this is the _sRGB_
    /// color space, which is rather small.
    pub(crate) space: RgbSpace,

    /// Reference to the colorimetric observer being used. This is almost always
    /// the CIE 1931 standard observer, which has been known to represent the
    /// deep blue region of humen vision sensitivity incorrectly. Here we allow
    /// other standard observers, such as the CIE 2015 cone fundamentals based
    /// observer, to improve color management quality.
    pub(crate) observer: Observer,
    pub(crate) rgb: Vector3<f64>,
}

impl RGB {
    /// Construct a RGB instance from red, green, and blue values in the range from 0 to 1,
    /// with an `Observer` and an optional `RgbSpace`.
    /// When using online RGB data, when observer and color space or color profile are not explicititely specfied,
    /// `Observer::Std1931`, and `RgbSpace::SRGB` are implied, and those are the defaults here too.
    pub fn new(
        r: f64,
        g: f64,
        b: f64,
        observer: Option<Observer>,
        space: Option<RgbSpace>,
    ) -> Self {
        let observer = observer.unwrap_or_default();
        let space = space.unwrap_or_default();
        RGB {
            rgb: Vector3::new(r, g, b),
            observer,
            space,
        }
    }

    /// Construct a RGB instance from red, green, and blue u8 values in the range from 0 to 255.
    ///
    /// When using online RGB data, when observer and color space or color profile are not explicititely specfied,
    /// `Observer::Std1931`, and `RgbSpace::SRGB` are implied, and those are the defaults here too.
    pub fn from_u8(
        r_u8: u8,
        g_u8: u8,
        b_u8: u8,
        observer: Option<Observer>,
        space: Option<RgbSpace>,
    ) -> Self {
        let space = space.unwrap_or_default();
        let [r, g, b] = [r_u8, g_u8, b_u8]
            .map(|v| (v as f64 / 255.0).clamp(0.0, 1.0))
            .map(|v| space.data().gamma.decode(v));
        RGB::new(r, g, b, observer, Some(space))
    }

    /// Construct a RGB instance from red, green, and blue u16 values in the range from 0 to 1.
    ///
    /// When using online RGB data, when observer and color space or color profile are not explicititely specfied,
    /// `Observer::Std1931`, and `RgbSpace::SRGB` are implied, and those are the defaults here too.
    pub fn from_u16(
        r_u16: u16,
        g_u16: u16,
        b_u16: u16,
        observer: Option<Observer>,
        space: Option<RgbSpace>,
    ) -> Self {
        let space = space.unwrap_or_default();
        let [r, g, b] = [r_u16, g_u16, b_u16]
            .map(|v| (v as f64 / 65_535.0).clamp(0.0, 1.0))
            .map(|v| space.data().gamma.decode(v));
        RGB::new(r, g, b, observer, Some(space))
    }

    pub fn from_xyz(xyz: XYZ, space: RgbSpace) -> Self {
        let xyz2rgb = xyz.observer.data().xyz2rgb(space);
        let xyz0 = xyz.xyz.unwrap_or(xyz.xyzn);
        let &[r, g, b] = (xyz2rgb * xyz0).as_ref();
        RGB::new(r, g, b, Some(xyz.observer), Some(space))
    }

    /// Returns the RGB values as an array with the red, green, and blue values respectively
    ///
    /// ```rust
    /// # use colorimetry::rgb::RGB;
    /// let rgb = RGB::new(0.1, 0.2, 0.3, None, None);
    /// let [r, g, b] = rgb.values();
    /// assert_eq!([r, g, b], [0.1, 0.2, 0.3]);
    /// ```
    pub fn values(&self) -> [f64; 3] {
        *self.rgb.as_ref()
    }

    /// Converts the RGB value to a tri-stimulus XYZ value
    pub fn xyz(&self) -> XYZ {
        const YW: f64 = 100.0;
        let xyzn = self
            .observer
            .data()
            .xyz(&self.space.data().white, None)
            .set_illuminance(100.0)
            .xyzn;
        let xyz = self.observer.data().rgb2xyz(&self.space) * self.rgb;
        XYZ {
            observer: self.observer,
            xyz: Some(xyz.map(|v| v * YW)),
            xyzn,
        }
    }

    // /// Creates a callback of closure function, which takes a set or RGB values, within a color
    // /// space and viewed as one observer, and returns a new set of RGB values, represeting the
    // /// stimulus in another color space, and using another observer.
    // ///
    // /// This conversion uses the spectral represenations of the primaries through the color space
    // /// `Spectra` function, to create a  transformation matrix.
    // pub fn convert(
    //     obs_from: Observer,
    //     space_from: RgbSpace,
    //     obs: Observer,
    //     space: RgbSpace,
    // ) -> Box<dyn Fn(&Vector3<f64>) -> Vector3<f64>> {
    //     todo!()
    // }

    // /// Transform a set of RGB values, defining a stimulus for one standard observer, into a set of
    // /// RGB values representing the same stimulus for different standard observer or special
    // /// observer.  On initial use this function calculates a transformation matrix based on the
    // /// colorimetric tristimulus values of the respective primaries.
    // pub fn transform(&self, obs_from: &Observer) -> Self {
    //     todo!()
    // }

    /*
    /// Creates a [Spectrum] from an [RGB] value, using the spectral primaries of its color space.
    /// See also [Spectrum::rgb] and [Spectrum::srgb].
    pub fn spectrum(&self) -> Spectrum {
        let p = &self.space.data().0.primaries;
        let yrgb = CIE1931.rgb2xyz(&RgbSpace::SRGB).row(1);
        self.data.iter().zip(yrgb.iter()).zip(p.iter()).map(|((v,w),s)|*v * *w * *s).sum::<Spectrum>().set_category(Category::Stimulus)
    }
     */
}

impl Light for RGB {
    fn spectrum(&self) -> Cow<Spectrum> {
        let prim = &self.space.data().primaries;
        let yrgb = self.observer.data().rgb2xyz(&self.space).row(1);
        //        self.rgb.iter().zip(yrgb.iter()).zip(prim.iter()).map(|((v,w),s)|*v * *w * &s.0).sum()
        let s = self
            .rgb
            .iter()
            .zip(yrgb.iter())
            .zip(prim.iter())
            .fold(Spectrum::default(), |acc, ((&v, &w), s)| acc + v * w * s.0);
        Cow::Owned(s)
    }
}

impl Filter for RGB {
    /**
        An RGB pixel as a filter.

        This excludes the reference white light.
        Ii is the filter function only, which is used in combination with a reference illuminant to achieve
        a stimulus in accordance with the colorspace in which is defined.

        ```
        use colorimetry::prelude::*;

        // rgb white in using CIE1931 standard observer, and sRGB color space.
        let rgb = RGB::from_u8(255, 255, 255, None, None);
        let d65: XYZ = CIE1931.xyz(&StdIlluminant::D65, Some(&rgb));
        approx::assert_ulps_eq!(d65, XYZ_D65WHITE, epsilon=1E-2);
        ```
    */
    fn spectrum(&self) -> Cow<Spectrum> {
        let prim = self.space.data().primaries_as_colorants();
        let yrgb = self.observer.data().rgb2xyz(&self.space).row(1);
        let s = self
            .rgb
            .iter()
            .zip(yrgb.iter())
            .zip(prim.iter())
            .fold(Spectrum::default(), |acc, ((&v, &w), s)| acc + v * w * s.0);
        Cow::Owned(s)
    }
}

impl AsRef<Vector3<f64>> for RGB {
    fn as_ref(&self) -> &Vector3<f64> {
        &self.rgb
    }
}

/// Clamped RGB values as a u8 array. Uses gamma function.
impl From<RGB> for [u8; 3] {
    fn from(rgb: RGB) -> Self {
        let data: &[f64; 3] = rgb.rgb.as_ref();
        data.map(|v| (rgb.space.data().gamma.encode(v.clamp(0.0, 1.0)) * 255.0).round() as u8)
    }
}

impl From<RGB> for [f64; 3] {
    fn from(rgb: RGB) -> Self {
        rgb.values()
    }
}

impl AbsDiffEq for RGB {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.observer == other.observer && self.rgb.abs_diff_eq(&other.rgb, epsilon)
    }
}

impl approx::UlpsEq for RGB {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.observer == other.observer && self.rgb.ulps_eq(&other.rgb, epsilon, max_ulps)
    }
}

pub fn gaussian_filtered_primaries(
    white: &Spectrum,
    red: [f64; 3],
    green: [f64; 2],
    blue: [f64; 2],
) -> [Stimulus; 3] {
    let [rc, rw, f] = red;
    let [gc, gw] = green;
    let [bc, bw] = blue;
    [
        Stimulus(
            Stimulus(&*Colorant::gaussian(bc, bw).spectrum() * white)
                .set_luminance(&CIE1931, 100.0)
                .0
                * f
                + Stimulus(&*Colorant::gaussian(rc, rw).spectrum() * white)
                    .set_luminance(&CIE1931, 100.0)
                    .0
                    * (1.0 - f),
        ),
        Stimulus(&*Colorant::gaussian(gc, gw).spectrum() * white).set_luminance(&CIE1931, 100.0),
        Stimulus(&*Colorant::gaussian(bc, bw).spectrum() * white).set_luminance(&CIE1931, 100.0),
    ]
}

#[cfg(test)]
mod rgb_tests {
    use crate::prelude::*;

    #[test]
    fn get_values_f64() {
        let rgb = RGB::new(0.1, 0.2, 0.3, None, None);
        let [r, g, b] = <[f64; 3]>::from(rgb);
        assert_eq!(r, 0.1);
        assert_eq!(g, 0.2);
        assert_eq!(b, 0.3);
    }

    #[test]
    fn get_values_u8() {
        let rgb = RGB::new(0.1, 0.2, 0.3, None, None);
        let [r, g, b] = <[u8; 3]>::from(rgb);
        assert_eq!(r, 89);
        assert_eq!(g, 124);
        assert_eq!(b, 149);
    }
}
