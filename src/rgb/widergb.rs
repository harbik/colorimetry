// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

//! # WideRgb: Color Representation allowing out-of-gamut colors for a given RGB color space.
//!
//! This module provides the `WideRgb` struct, a representation of a color stimulus using unconstrained
//! Red, Green, and Blue (RGB) floating-point values within a specified RGB color space. Unlike typical
//! RGB values that are limited to the range `[0.0, 1.0]`, `WideRgb` allows values to extend beyond this
//! range, accommodating out-of-gamut colors that cannot be accurately rendered by most devices.
//!
//! ## Key Features
//! - **Unconstrained RGB Values:** The `WideRgb` struct stores RGB values that can extend beyond the
//!   standard range, allowing for high-dynamic-range (HDR) and out-of-gamut colors.
//! - **Color Space Support:** Supports various RGB color spaces, including sRGB, Adobe RGB, and custom
//!   spaces, specified through the `RgbSpace` type.
//! - **Observer Customization:** Allows the use of different colorimetric observers (e.g., CIE 1931,
//!   CIE 2015), enhancing accuracy in color conversions and comparisons.
//! - **Conversion Methods:** Provides methods to clamp, compress, and convert `WideRgb` values to
//!   standard `Rgb` values within the `[0.0, 1.0]` range.
//! - **Compatibility:** Includes conversions to `XYZ` tristimulus values and standard `Rgb`.
//!
//! ## Example Usage
//!
//! ```rust
//! use colorimetry::rgb::WideRgb;
//!
//! let wide_rgb = WideRgb::new(1.5, -0.2, 0.8, None, None);
//! let compressed_rgb = wide_rgb.compress();
//! let clamped_rgb = wide_rgb.clamp();
//!
//! println!("Compressed: {:?}", compressed_rgb.to_array());
//! println!("Clamped: {:?}", clamped_rgb.to_array());
//! ```
//!
//! ## Notes
//! - The `compress` method applies chroma reduction and luminance scaling to ensure values fit within
//!   the `[0.0, 1.0]` range, resulting in a lossy conversion.
//! - The `clamp` method simply restricts values to the valid range without chroma adjustment, potentially
//!   altering perceived color balance.
//!
//! ## Testing
//! Unit tests are provided to verify the accuracy of conversions, clamping, and compression operations.

use crate::{
    observer::Observer,
    rgb::{rgbspace::RgbSpace, Rgb},
    xyz::XYZ,
};
use approx::AbsDiffEq;
use nalgebra::Vector3;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Debug, Clone, Copy, PartialEq)]
/// Represents a color stimulus using unconstrained Red, Green, and Blue (RGB) floating-point values
/// within a device's RGB color space. The values can extend beyond the typical 0.0 to 1.0 range,
/// allowing for out-of-gamut colors that cannot be accurately represented by the device.
pub struct WideRgb {
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

impl WideRgb {
    /// Creates a new `WideRgb` instance with the specified red, green, and blue values.
    /// # Parameters
    /// - `r`: Red component, any f64 value
    /// - `g`: Green component, any f64 value
    /// - `b`: Blue component, any f64 value
    /// - `observer`: Optional observer, defaults to `Observer::Cie1931`
    /// - `space`: Optional RGB color space, defaults to `RgbSpace::SRGB`
    /// # Returns
    /// A new `WideRgb` instance with the specified RGB values and color space.
    pub fn new(
        r: f64,
        g: f64,
        b: f64,
        observer: Option<Observer>,
        space: Option<RgbSpace>,
    ) -> Self {
        let observer = observer.unwrap_or_default();
        let space = space.unwrap_or_default();
        WideRgb {
            rgb: Vector3::new(r, g, b),
            observer,
            space,
        }
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
    /// # use colorimetry::rgb::WideRgb;
    /// let rgb = WideRgb::new(0.1, 0.2, 0.3, None, None);
    /// let [r, g, b] = rgb.to_array();
    /// assert_eq!([r, g, b], [0.1, 0.2, 0.3]);
    /// ```
    pub fn to_array(&self) -> [f64; 3] {
        *self.rgb.as_ref()
    }

    /// Converts the RGB value to a tri-stimulus XYZ value
    ///
    /// This method applies the observer's RGB to XYZ conversion matrix to the RGB values,
    /// resulting in XYZ values scaled by the white point's Y value (100.0).
    ///  # Returns
    /// An `XYZ` instance representing the color in the XYZ color space.
    /// ```
    /// # use colorimetry::rgb::WideRgb;
    /// # use approx::assert_abs_diff_eq;
    /// # use colorimetry::xyz::XYZ;
    /// # use colorimetry::observer::Observer::Cie1931;
    /// # use colorimetry::rgb::RgbSpace::SRGB;
    /// let xyz = Cie1931.xyz_d65();
    /// let wide_rgb = xyz.rgb(SRGB);
    /// let xyz_back = wide_rgb.xyz();
    /// assert_abs_diff_eq!(xyz, xyz_back, epsilon = 5E-4); // Example assertion, actual values depend on the observer and space
    /// ```
    ///
    pub fn xyz(&self) -> XYZ {
        const YW: f64 = 100.0;
        let xyz = self.observer.rgb2xyz_matrix(self.space) * self.rgb;
        XYZ {
            observer: self.observer,
            xyz: xyz.map(|v| v * YW),
        }
    }

    /// Returns whether or not this Wide RGB instance is within the RGB gamut or not.
    ///
    /// Being in gamut means all color channel values are in the range [0.0, 1.0].
    ///
    /// ```rust
    /// # use colorimetry::rgb::WideRgb;
    ///
    /// let in_gamut = WideRgb::new(1.0, 0.2, 0.8, None, None);
    /// assert!(in_gamut.is_in_gamut());
    ///
    /// let out_of_gamut = WideRgb::new(1.2, -0.5, 0.8, None, None);
    /// assert!(!out_of_gamut.is_in_gamut());
    /// ```
    pub fn is_in_gamut(&self) -> bool {
        self.to_array().iter().all(|v| (0.0..=1.0).contains(v))
    }

    /// Returns whether or not this Wide RGB instance is within the RGB gamut or not,
    /// excluding the border values. Use `epsilon` to define the margin
    /// ```rust
    /// # use colorimetry::rgb::WideRgb;
    ///
    /// let in_gamut = WideRgb::new(0.0, 0.3, 0.3, None, None);
    /// assert!(in_gamut.is_within_gamut(0.001));
    ///
    /// let in_gamut = WideRgb::new(0.99, 0.3, 0.3, None, None);
    /// assert!(in_gamut.is_within_gamut(0.001));
    ///
    /// let in_gamut = WideRgb::new(0.5, 0.5, 0.9999, None, None);
    /// assert!(!in_gamut.is_within_gamut(0.01));
    ///
    /// let in_gamut = WideRgb::new(0.5, 0.5, 0.989, None, None);
    /// assert!(in_gamut.is_within_gamut(0.01));
    /// ```
    pub fn is_within_gamut(&self, epsilon: f64) -> bool {
        self.to_array().iter().all(|&v| v < 1.0 - epsilon)
    }

    pub fn is_black(&self, epsilon: f64) -> bool {
        self.rgb.iter().all(|&v| v.abs() < epsilon)
    }

    /// Returns itself as an [`Rgb`] instance with all channel values unchanged, if this wide RGB
    /// instance is in the RGB gamut (all channel values within the `[0.0, 1.0]` range).
    /// If this wide RGB instance is out-of-gamut, `None` is returned.
    ///
    /// ```rust
    /// # use colorimetry::rgb::WideRgb;
    ///
    /// let wide_rgb = WideRgb::new(1.0, 0.2, 0.8, None, None);
    /// let rgb = wide_rgb.to_rgb().unwrap();
    /// assert_eq!(rgb.to_array(), wide_rgb.to_array());
    ///
    /// // A `WideRgb` value with out-of-gamut components.
    /// let wide_rgb = WideRgb::new(1.2, -0.5, 0.8, None, None);
    /// assert!(wide_rgb.to_rgb().is_none());
    /// ```
    pub fn to_rgb(self) -> Option<Rgb> {
        if self.is_in_gamut() {
            Some(Rgb {
                rgb: self.rgb,
                observer: self.observer,
                space: self.space,
            })
        } else {
            None
        }
    }

    /// Converts a `WideRgb` value to a valid `Rgb` value by clamping red, green, and blue values to the range [0, 1].
    ///
    /// ```rust
    /// # use colorimetry::rgb::WideRgb;
    ///
    /// // A `WideRgb` value with out-of-gamut components.
    /// let wide_rgb = WideRgb::new(1.2, -0.5, 0.8, None, None);
    ///
    /// // Clamp the values to the [0.0, 1.0] range.
    /// let clamped_rgb = wide_rgb.clamp();
    ///
    /// assert_eq!(clamped_rgb.to_array(), [1.0, 0.0, 0.8]);
    /// ```
    ///
    /// # Parameters
    /// - `self`: The `WideRgb` instance to be clamped.
    ///
    /// # Returns
    /// A new `Rgb` instance with the adjusted RGB values, ensuring they are within the allowable range.
    pub fn clamp(&self) -> Rgb {
        let rgb = self.rgb.map(|v| v.clamp(0.0, 1.0));

        Rgb {
            rgb,
            observer: self.observer,
            space: self.space,
        }
    }

    /// Converts a `WideRgb` value to a valid `Rgb` value by compressing out-of-gamut colors, if necessary.
    /// This method adjusts out-of-gamut colors to fit within the RGB color space by reducing chroma (adding white)
    /// and scaling luminance to the device’s maximum allowable value.
    ///
    /// # Parameters
    /// - `self`: The `WideRgb` instance to be processed. If the values are already within the RGB gamut,
    ///   they will be returned unchanged. Otherwise, the values will be compressed to fit within the valid range.
    ///
    /// # Returns
    /// A new `Rgb` instance with out-of-gamut colors compressed to fit within the RGB color space, ensuring
    /// that all colors can be accurately rendered by the device.
    ///
    /// # Example
    /// ```rust
    /// # use colorimetry::rgb::WideRgb;
    /// # use colorimetry::rgb::Rgb;
    /// # use approx::assert_abs_diff_eq;
    ///
    /// // A `WideRgb` value with out-of-gamut components.
    /// let wide_rgb = WideRgb::new(1.2, -0.5, 0.8, None, None);
    ///
    /// // Compresses the values to the [0.0, 1.0] range.
    /// let compressed_rgb = wide_rgb.compress();
    ///
    /// assert_abs_diff_eq!(compressed_rgb.to_array().as_ref(), [1.0, 0.0, 0.7647].as_ref(), epsilon = 0.0001);
    /// ```
    ///
    /// # Notes
    /// This is a lossy operation. Out-of-gamut colors are modified by changing chroma and luminance values,
    /// resulting in a color that differs from the original.
    pub fn compress(&self) -> Rgb {
        // Amount to add to get all channels positive
        let translate = -self.rgb.min().min(0.0);
        // The scaling needed to get all channels below 1.0
        let scale = (self.rgb.max() + translate).max(1.0);

        let in_gamut_rgb = if translate != 0.0 || scale != 1.0 {
            self.rgb.add_scalar(translate) / scale
        } else {
            self.rgb
        };
        Rgb {
            rgb: in_gamut_rgb,
            observer: self.observer,
            space: self.space,
        }
    }
}

impl AsRef<Vector3<f64>> for WideRgb {
    fn as_ref(&self) -> &Vector3<f64> {
        &self.rgb
    }
}

impl From<WideRgb> for [f64; 3] {
    fn from(rgb: WideRgb) -> Self {
        rgb.to_array()
    }
}

impl AbsDiffEq for WideRgb {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.observer == other.observer && self.rgb.abs_diff_eq(&other.rgb, epsilon)
    }
}

impl approx::UlpsEq for WideRgb {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.observer == other.observer && self.rgb.ulps_eq(&other.rgb, epsilon, max_ulps)
    }
}

#[cfg(test)]
mod rgb_tests {
    use crate::rgb::WideRgb;

    #[test]
    fn get_values() {
        let rgb = WideRgb::new(0.1, 0.2, 0.3, None, None);
        let [r, g, b] = rgb.to_array();
        assert_eq!(r, 0.1);
        assert_eq!(g, 0.2);
        assert_eq!(b, 0.3);

        assert_eq!(rgb.r(), 0.1);
        assert_eq!(rgb.g(), 0.2);
        assert_eq!(rgb.b(), 0.3);

        assert_eq!(<[f64; 3]>::from(rgb), rgb.to_array());
    }

    #[test]
    fn out_of_gamut_values() {
        let rgb = WideRgb::new(-0.8, 2.7, 0.8, None, None);
        let [r, g, b] = rgb.to_array();
        assert_eq!(r, -0.8);
        assert_eq!(g, 2.7);
        assert_eq!(b, 0.8);
    }

    #[test]
    fn compress_in_gamut() {
        let rgb = WideRgb::new(0.0, 0.0, 0.0, None, None).compress();
        assert_eq!(rgb.to_array(), [0.0, 0.0, 0.0]);

        let rgb = WideRgb::new(1.0, 1.0, 1.0, None, None).compress();
        assert_eq!(rgb.to_array(), [1.0, 1.0, 1.0]);

        let rgb = WideRgb::new(0.75, 0.75, 0.75, None, None).compress();
        assert_eq!(rgb.to_array(), [0.75, 0.75, 0.75]);

        let rgb = WideRgb::new(0.1, 0.3, 1.0, None, None).compress();
        assert_eq!(rgb.to_array(), [0.1, 0.3, 1.0]);

        let rgb = WideRgb::new(0.0, 0.3, 1.0, None, None).compress();
        assert_eq!(rgb.to_array(), [0.0, 0.3, 1.0]);
    }

    #[test]
    fn compress_out_of_gamut() {
        // Same value above 1.0 are all scaled down to 1.0
        let rgb = WideRgb::new(1.2, 1.2, 1.2, None, None).compress();
        assert_eq!(rgb.to_array(), [1.0, 1.0, 1.0]);

        // Same value below 0.0 are all scaled up to 0.0
        let rgb = WideRgb::new(-0.5, -0.5, -0.5, None, None).compress();
        assert_eq!(rgb.to_array(), [0.0, 0.0, 0.0]);

        // Different positiv values are all scaled down with the max channel value
        let rgb = WideRgb::new(1.0, 2.0, 3.0, None, None).compress();
        assert_eq!(rgb.to_array(), [1.0 / 3.0, 2.0 / 3.0, 1.0]);

        // Values outside the range both above and below are compressed correctly
        let rgb = WideRgb::new(-0.1, 0.9, 1.1, None, None).compress();
        assert_eq!(rgb.to_array(), [0.0, (0.9 + 0.1) / (1.1 + 0.1), 1.0]);

        // Some negative channel, and some channel that becomes too large after
        // adding the lowest channel to it.
        let rgb = WideRgb::new(-0.5, 0.4, 0.6, None, None).compress();
        assert_eq!(rgb.to_array(), [0.0, (0.4 + 0.5) / (0.6 + 0.5), 1.0]);

        // One negative channel, and the rest stays within gamut even after adding
        // the lowest channel to it.
        let rgb = WideRgb::new(-0.5, 0.4, 0.4, None, None).compress();
        assert_eq!(rgb.to_array(), [0.0, 0.9, 0.9]);
    }
}
