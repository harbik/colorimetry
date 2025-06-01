//! # Device-dependent Red-Green-Blue (RGB) Colors
//!
//! The `rgb` module provides types and functions for working with device-dependent
//! Red-Green-Blue colors.
//!
//! ## Submodules
//!
//! - **`gamma`**  
//!   Encoding and decoding gamma curves (e.g. sRGB, AdobeRGB) for nonlinear ↔ linear conversions.  
//! - **`rgbspace`**  
//!   Definitions of `RgbSpace` (primaries, white point, conversion matrices) and helpers  
//!   to transform between RGB and CIE XYZ.  
//! - **`widergb`**  
//!   The `WideRgb` type, which allows out-of-gamut values (outside `[0.0,1.0]`) for HDR or extended-gamut workflows.
//!
//! ## Core Features
//!
//! - **Strict Validation**  
//!   Creating an `Rgb` with `Rgb::new(r, g, b, …)` returns an error if any component lies
//!   outside `[0.0, 1.0]`.  
//! - **Multiple Constructors**  
//!   - `Rgb::new(r, g, b, observer, space)` — fully specified  
//!   - `Rgb::from_u8(r, g, b, …)` / `Rgb::from_u16(r, g, b, …)` — byte- or word-based input  
//! - **Color Space & Observer**  
//!   Optionally supply a standard observer (e.g. CIE1931 or CIE2015) and an RGB color
//!   space (e.g. sRGB, AdobeRGB). Defaults are Std1931 + sRGB.
//!
//! ## Conversions
//!
//! - **To XYZ**  
//!   Call `.xyz()` to convert into CIE XYZ tristimulus values (scaled so Y = 100 for white).
//! - **To Device Values**  
//!   Use `Into<[u8; 3]>` or `Into<[f64; 3]>` for byte or float arrays, automatically applying
//!   the color space’s gamma curve.
//!
//! ## Error Handling
//!
//! - `Rgb::new(...)` returns `Err(CmtError::InvalidRgbValue)` if any component ∉ `[0.0,1.0]`.  
//! - `from_u8` / `from_u16` automatically clamp inputs into range and never error.
//!
//! ## Notes
//!
//! - Use `Rgb` when you need **strict gamut compliance**. For HDR or out-of-gamut workflows,
//!   consider `WideRgb`.  
//! - Supplying an alternate `observer` or `space` customizes the `.xyz()` conversion for specialized
//!   viewing conditions or display profiles.
//!
//! ## Testing
//!
//! The module includes unit tests for:
//! - Value validation and clamping  
//! - Byte and float conversions  
//! - XYZ tristimulus output under default settings  

mod gamma;
mod rgbspace;
mod widergb;

pub use rgbspace::RgbSpace;
pub use widergb::WideRgb;

use crate::{
    colorant::Colorant,
    error::Error,
    observer::Observer,
    observer::CIE1931,
    spectrum::Spectrum,
    stimulus::Stimulus,
    traits::{Filter, Light},
    xyz::XYZ,
};
use approx::AbsDiffEq;
use nalgebra::Vector3;
use std::borrow::Cow;

/// Represents a color stimulus using Red, Green, and Blue (RGB) values constrained to the `[0.0, 1.0]` range.
/// Each component is a floating-point value representing the relative intensity of the respective primary color
/// within a defined RGB color space.
///
/// Unlike the CIE XYZ tristimulus values, which use imaginary primaries, RGB values are defined using real primaries
/// based on a specific color space. These primaries typically form a triangular area within a CIE (x,y) chromaticity
/// diagram, representing the gamut of colors the device can reproduce.
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
#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
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
    ) -> Result<Self, Error> {
        if (0.0..=1.0).contains(&r) && (0.0..=1.0).contains(&g) && (0.0..=1.0).contains(&b) {
            let observer = opt_observer.unwrap_or_default();
            let space = opt_rgbspace.unwrap_or_default();
            Ok(Rgb {
                rgb: Vector3::new(r, g, b),
                observer,
                space,
            })
        } else {
            Err(Error::InvalidRgbValue)
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

    /// Returns the value of the red channel.
    pub fn r(&self) -> f64 {
        self.rgb.x
    }

    /// Returns the value of the green channel.
    pub fn g(&self) -> f64 {
        self.rgb.y
    }

    /// Returns the value of the blue channel.
    pub fn b(&self) -> f64 {
        self.rgb.z
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
        let xyz = self.observer.data().rgb2xyz(&self.space) * self.rgb;
        XYZ {
            observer: self.observer,
            xyz: xyz.map(|v| v * YW),
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
    /// let d65: XYZ = CIE1931.xyz(&CieIlluminant::D65, Some(&rgb));
    /// let xyz_d65 = CIE1931.xyz_d65();
    /// approx::assert_ulps_eq!(d65, xyz_d65, epsilon = 1e-2);
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
