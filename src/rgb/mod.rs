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
//! use colorimetry::rgb::rgb::Rgb;
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

pub mod rgb;
pub mod rgbspace;
pub mod widergb;
