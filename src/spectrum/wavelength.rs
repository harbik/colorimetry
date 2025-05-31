//! # Functions and Constants from Physics
//!
//! This module provides core physical constants and spectral‐generation functions
//! used throughout the library to build illuminants, LEDs, black‐body radiators,
//! and band-pass filters.  
//!
//! ## Wavelength Conversions
//! - `wavelength(i)` — Map an index or meter‐value `i` to meters (handles both raw nm and index → nm).  
//! - `to_wavelength(x, xmin, xmax)` — Linearly interpolate a domain value `x ∈ [xmin, xmax]` to the
//!   spectrum range [380 nm … 780 nm].  
//!
//! ## Usage
//! Use these primitives to generate accurate spectral power distributions for illuminants,
//! stimuli, colorants, and to perform photometric and radiometric calculations.
//!

use crate::spectrum::SPECTRUM_WAVELENGTH_RANGE;
use core::f64;
use num_traits::ToPrimitive;

#[inline(always)]
pub fn wavelength<T: ToPrimitive>(i: T) -> f64 {
    let f = i.to_f64().unwrap_or(f64::NAN);
    if f > 1E-3 {
        f * 1E-9
    } else {
        f
    }
}

/// Map a value x, in a domain from xmin to xmax to a wavelength in the domain
/// from 380E-9 to 780E-9 meter.
/// ```
/// // Wavelength from an index value in the domain from 0 to 400:
/// use colorimetry::spectrum::to_wavelength;
///
/// let l = to_wavelength(200, 0, 400);
/// approx::assert_ulps_eq!(l, 580E-9);
///
/// // Wavelength from from a function defined over a domain from 0.0 to 1.0:
/// let l = to_wavelength(0.0, 0.0, 1.0);
/// approx::assert_ulps_eq!(l, 380E-9);
///
/// // Wavelength defined in integer nanometer values, to floating point meters
/// let l = to_wavelength(780, 380, 780);
/// approx::assert_ulps_eq!(l, 780E-9);
/// ```
#[inline]
pub fn to_wavelength<T: ToPrimitive>(x: T, xmin: T, xmax: T) -> f64 {
    let xmin = xmin.to_f64().unwrap_or(f64::NAN);
    let xmax = xmax.to_f64().unwrap_or(f64::NAN);
    let x = x.to_f64().unwrap_or(f64::NAN);
    let f = (x - xmin) / (xmax - xmin);
    wavelength(*SPECTRUM_WAVELENGTH_RANGE.start()) * (1.0 - f)
        + wavelength(*SPECTRUM_WAVELENGTH_RANGE.end()) * f
}

/// Convenience function for specifying wavelengths in nanometers or meters.
///
/// This accepts integer and float values.
/// Wwavelength values larger than 1E-3 are assumed to have the unit nanometer
/// and are converted to a unit of meters.
/// All integer values are nanometaer values.
pub fn wavelengths<T: ToPrimitive, const N: usize>(v: [T; N]) -> [f64; N] {
    v.map(|x| wavelength(x))
}
