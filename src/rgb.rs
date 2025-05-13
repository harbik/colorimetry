//! # Rgb: In gamut RGB Color Representation of a Color Stimulus
//!
//! This module provides the `Rgb` struct, a representation of a color stimulus using Red, Green, and Blue
//! (RGB) floating-point values constrained to the `[0.0, 1.0]` range. Unlike the `WideRgb` type, which allows
//! for out-of-gamut colors, `Rgb` is strictly limited to the device’s RGB color space gamut, ensuring that
//! all values are within the valid range for rendering and display.
//!
//! ## Key Features
//! - **Constrained RGB Values:** The `Rgb` struct enforces that all RGB values must be within the range
//!   `[0.0, 1.0]`. Attempts to create values outside this range will result in an error.
//! - **Color Space Support:** Supports various RGB color spaces through the `RgbSpace` type, including sRGB,
//!   Adobe RGB, and custom-defined spaces.
//! - **Observer Customization:** Allows the use of different colorimetric observers (e.g., CIE 1931, CIE 2015),
//!   enhancing accuracy in color conversions.
//! - **Data Conversions:** Provides methods to convert RGB values to `XYZ` tristimulus values, and to u8 and u16
//!   arrays for compatibility with 8-bit and 16-bit RGB formats.
//! - **Validation:** Ensures that all RGB values are validated during creation, preventing invalid color data.
//!
//! ## Example Usage
//!
//! ```rust
//! use colorimetry::rgb::Rgb;
//!
//! // Creating a valid Rgb instance
//! let rgb = Rgb::new(0.5, 0.3, 0.7, None, None).unwrap();
//!
//! // Convert to a u8 array for image processing
//! let u8_rgb: [u8; 3] = rgb.into();
//!
//! // Convert to XYZ tristimulus values
//! let xyz = rgb.xyz();
//!
//! println!("RGB as u8: {:?}", u8_rgb);
//! println!("XYZ values: {:?}", xyz.values());
//! ```
//!
//! ## Error Handling
//! - The `Rgb::new()` method returns an `Err(CmtError::InvalidRgbValue)` if any of the provided RGB values
//!   are outside the `[0.0, 1.0]` range.
//! - The `from_u8()` and `from_u16()` methods internally clamp values to ensure they remain within the valid range.
//!
//! ## Notes
//! - Unlike the `WideRgb` type, which supports out-of-gamut colors, `Rgb` enforces strict adherence to the
//!   `[0.0, 1.0]` range. As such, it is ideal for cases where data must be strictly constrained to the device’s
//!   rendering gamut.
//! - The `observer` field is optional but can be specified to provide more accurate conversions to `XYZ` values
//!   based on specific viewing conditions.
//!
//! ## Testing
//! - Comprehensive unit tests verify RGB value validation, conversion methods, and compatibility with different
//!   data formats, including `u8` and `f64` arrays.

use crate::{
    colorant::Colorant,
    data::observers::CIE1931,
    error::CmtError,
    illuminant::Illuminant,
    observer::Observer,
    prelude::WideRgb,
    rgbspace::RgbSpace,
    spectrum::Spectrum,
    stimulus::Stimulus,
    traits::{Filter, Light},
    xyz::XYZ,
};
use approx::AbsDiffEq;
use nalgebra::{Matrix3, Vector3};
use std::borrow::Cow;
use wasm_bindgen::prelude::wasm_bindgen;

/// Represents a color stimulus using Red, Green, and Blue (RGB) values constrained to the `[0.0, 1.0]` range.
/// Each component is a floating-point value representing the relative intensity of the respective primary color
/// within a defined RGB color space.
///
/// Unlike the CIE XYZ tristimulus values, which use imaginary primaries, RGB values are defined using real primaries
/// based on a specific color space. These primaries typically form a triangular area within a CIE (x,y) chromaticity
/// diagram, representing the gamut of colors the device can reproduce.
///
/// # Fields
/// - `space`: The RGB color space defining the specific primaries and white point. Common spaces include sRGB, Adobe RGB, and custom-defined spaces.
/// - `observer`: The colorimetric observer that defines how the color is perceived. Default is the CIE 1931 standard observer,
///   but more recent observers, such as the CIE 2015 cone fundamentals, can be specified for improved accuracy.
/// - `rgb`: A 3-element vector representing the red, green, and blue components as floating-point values in the range `[0.0, 1.0]`.
///
/// # Usage
/// The `Rgb` struct is used to encapsulate color information in a device-independent manner, allowing for accurate color
/// representation, conversion, and manipulation within defined RGB spaces. It is particularly useful for applications
/// involving color management, digital imaging, and rendering where strict adherence to gamut boundaries is required.
///
/// # Example
/// ```rust
/// # use colorimetry::rgb::Rgb;
/// # use approx::assert_abs_diff_eq;
///
/// // Create an sRGB color with normalized RGB values
/// let rgb = Rgb::new(0.5, 0.25, 0.75, None, None).unwrap();
/// assert_abs_diff_eq!(rgb.values().as_ref(), [0.5, 0.25, 0.75].as_ref(), epsilon = 1e-6);
/// ```
///
/// # Notes
/// - The `Rgb` struct strictly enforces the `[0.0, 1.0]` range for each component. Any attempt to create values
///   outside this range will result in an error.
/// - The `observer` field allows for color conversion accuracy under different lighting and viewing conditions,
///   enhancing the reliability of transformations to other color spaces such as XYZ.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rgb {
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

impl Rgb {
    /// Creates a new `Rgb` instance with the specified red, green, and blue values.
    /// # Arguments
    /// - `r`: Red component, in the range from 0.0 to 1.0
    /// - `g`: Green component, in the range from 0.0 to 1.0
    /// - `b`: Blue component, in the range from 0.0 to 1.0
    /// - `observer`: Optional observer, defaults to `Observer::Std1931`
    /// - `space`: Optional RGB color space, defaults to `RgbSpace::SRGB`
    /// # Returns
    /// A new `Rgb` instance with the specified RGB values and color space.
    /// # Errors
    /// - InvalidRgbValue: at least one of the RGB values is outside the range from 0.0 to 1.0.
    ///     
    pub fn new(
        r: f64,
        g: f64,
        b: f64,
        opt_observer: Option<Observer>,
        opt_rgbspace: Option<RgbSpace>,
    ) -> Result<Self, CmtError> {
        if r >= 0.0 && r <= 1.0 && g >= 0.0 && g <= 1.0 && b >= 0.0 && b <= 1.0 {
            let observer = opt_observer.unwrap_or_default();
            let space = opt_rgbspace.unwrap_or_default();
            Ok(Rgb {
                rgb: Vector3::new(r, g, b),
                observer,
                space,
            })
        } else {
            Err(CmtError::InvalidRgbValue)
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
        // unwrap OK as derived from u8 values and clamped to 0.0..1.0
        Rgb::new(r, g, b, observer, Some(space)).unwrap()
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
        // unwrap OK as derived from u16 values and clamped to 0.0..1.0
        Rgb::new(r, g, b, observer, Some(space)).unwrap()
    }

    /// Returns the RGB values as an array with the red, green, and blue values respectively
    ///
    /// ```rust
    /// # use colorimetry::rgb::Rgb;
    /// let rgb = Rgb::new(0.1, 0.2, 0.3, None, None).unwrap();
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
}


impl Light for Rgb {

    /// Implements the `Light` trait for the `Rgb` struct, allowing an `Rgb` color to be interpreted
    /// as a spectral power distribution (`Spectrum`).
    ///
    /// The `spectrum` method converts the RGB values into a spectral representation based on the
    /// primaries defined by the associated RGB color space and the current observer.
    ///
    /// # Method: `spectrum()`
    /// - This method calculates the spectral power distribution of the `Rgb` instance by combining the
    ///   RGB values with the respective primary spectra.
    /// - Each RGB component is weighted by its corresponding luminance factor (`yrgb`) and primary spectrum,
    ///   effectively converting the RGB values into a continuous spectrum representation.
    ///
    /// # Implementation Details
    /// - The method iterates over the RGB components and the respective primary spectra.
    /// - Each component value is scaled by its corresponding luminance weight (`yrgb`) and then combined
    ///   with the primary spectrum using a weighted sum.
    /// - The resulting `Spectrum` is returned as an owned `Cow<Spectrum>`.
    ///
    /// # Notes
    /// - The spectral representation is device-dependent and based on the primaries defined by the `RgbSpace`.
    /// - The observer's data is used to apply luminance scaling, enhancing perceptual accuracy.
    fn spectrum(&self) -> Cow<Spectrum> {
        let prim = &self.space.data().primaries;
        let rgb2xyz = self.observer.data().rgb2xyz(&self.space);
        let yrgb = rgb2xyz.row(1);
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

impl Filter for Rgb {
    /// Implements the `Filter` trait for the `Rgb` struct, treating an RGB color as a spectral filter function.
    ///
    /// The `spectrum` method interprets the RGB values as a filter, excluding the reference illuminant. 
    /// The filter function is then used in combination with a reference illuminant to simulate the resulting 
    /// stimulus within the defined RGB color space.
    ///
    /// # Method: `spectrum()`
    /// - This method calculates the spectral power distribution of the `Rgb` instance by treating each RGB 
    ///   component as a filter over its respective primary spectrum.
    /// - The resulting spectrum represents the relative transmittance of each primary, scaled by its respective 
    ///   luminance weight (`yrgb`), without applying any reference illuminant.
    ///
    /// # Example
    /// ```rust
    /// use colorimetry::prelude::*;
    /// use approx::assert_ulps_eq;
    /// 
    /// // Define an sRGB white color using the CIE 1931 observer
    /// let rgb = Rgb::from_u8(255, 255, 255, None, None);
    ///
    /// // Retrieve the spectral representation
    /// let spectrum = Filter::spectrum(&rgb);
    ///
    /// // Compare with the CIE D65 reference white point
    /// let d65: XYZ = CIE1931.xyz(&StdIlluminant::D65, Some(&rgb));
    /// approx::assert_ulps_eq!(d65, XYZ_D65WHITE, epsilon = 1e-2);
    /// ```
    ///
    /// # Implementation Details
    /// - The method iterates over the RGB components and their respective primary spectra, treating each 
    ///   component as a filter function.
    /// - Each component is scaled by its luminance factor (`yrgb`) to accurately reflect the relative 
    ///   contribution of each primary to the resulting spectrum.
    /// - The resulting spectrum is returned as an owned `Cow<Spectrum>`.
    ///
    /// # Notes
    /// - The spectral representation is device-dependent, relying on the primary spectra defined in the 
    ///   associated `RgbSpace`.
    /// - This implementation excludes the reference illuminant, making it suitable for use as a relative filter 
    ///   that can be combined with any illuminant to produce a specific stimulus.
    fn spectrum(&self) -> Cow<Spectrum> {
        let prim = self.space.data().primaries_as_colorants();
        let rgb2xyz = self.observer.data().rgb2xyz(&self.space);
        let yrgb = rgb2xyz.row(1);
        let s = self
            .rgb
            .iter()
            .zip(yrgb.iter())
            .zip(prim.iter())
            .fold(Spectrum::default(), |acc, ((&v, &w), s)| acc + v * w * s.0);
        Cow::Owned(s)
    }
}

impl AsRef<Vector3<f64>> for Rgb {
    fn as_ref(&self) -> &Vector3<f64> {
        &self.rgb
    }
}

/// Clamped RGB values as a u8 array. Uses gamma function.
impl From<Rgb> for [u8; 3] {
    fn from(rgb: Rgb) -> Self {
        let data: &[f64; 3] = rgb.rgb.as_ref();
        data.map(|v| (rgb.space.data().gamma.encode(v.clamp(0.0, 1.0)) * 255.0).round() as u8)
    }
}

impl From<Rgb> for [f64; 3] {
    fn from(rgb: Rgb) -> Self {
        rgb.values()
    }
}

impl AbsDiffEq for Rgb {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.observer == other.observer && self.rgb.abs_diff_eq(&other.rgb, epsilon)
    }
}

impl approx::UlpsEq for Rgb {
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
        let rgb = Rgb::new(0.1, 0.2, 0.3, None, None).unwrap();
        let [r, g, b] = <[f64; 3]>::from(rgb);
        assert_eq!(r, 0.1);
        assert_eq!(g, 0.2);
        assert_eq!(b, 0.3);
    }

    #[test]
    fn get_values_u8() {
        let rgb = Rgb::new(0.1, 0.2, 0.3, None, None).unwrap();
        let [r, g, b] = <[u8; 3]>::from(rgb);
        assert_eq!(r, 89);
        assert_eq!(g, 124);
        assert_eq!(b, 149);
    }
}
