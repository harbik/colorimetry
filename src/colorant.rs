//! # Spectral Reflectance/Transmission Functions for Filters and Surface Colors
//!
//! The **Colorant** module defines spectral **reflectance/transmission filters** (colorants or color patches),
//! represented as `Spectrum` values in the range [0.0…1.0] over 380 nm…780 nm at 1 nm steps (401 samples).
//! A `Colorant` can model physical filters, paints, or ideal‐theoretical patches.
//!
//! ## Key Type
//!
//! - `pub struct Colorant`
//!   Wraps a `Spectrum` of 401 floating-point reflectance values, clamped to \[0,1\].  
//!
//! ## Constructors & Factories
//!
//! ```text
//! Colorant::new(spectrum)        // Validate and wrap a raw Spectrum
//! Colorant::gray(g)              // Uniform reflectance = g
//! Colorant::white()              // Perfect white patch (guu=1.0)
//! Colorant::black()              // Perfect black patch (g=0.0)
//! Colorant::top_hat(center, w)   // Rectangular bandpass filter
//! Colorant::gaussian(center, σ)  // Gaussian band filter
//! ```
//!
//! You can also construct analytically:
//! ```text
//! let f: Colorant = (|x: f64| /* f(x) in [0,1] */).into();
//! ```
//!
//! ## Conversions & Metrics
//!
//! - **CIELab**:  
//!   ```text
//!   let lab: CieLab = colorant.cielab(Some(&illuminant), Some(observer));
//!   ```  
//!   Converts reflectance+illuminant → tristimulus → L\*a\*b\*.
//!
//! ## Arithmetic & Traits
//!
//! - **Subtractive mixing** via `*` (multiplicative spectrum clamp)  
//! - **Additive combination** via `+` (spectrum clamp)  
//! - Implements `Filter` and `Light` for seamless use in `Observer::xyz`  
//! - Supports `Mul<f64>`, `AddAssign`, `AbsDiffEq` for testing and scaling  
//!
//! ## Optional Features
//!
//! - `#[cfg(feature = "munsell")]`  
//!   Includes the `munsell_matt` module for Munsell‐based matte patches.
//!

#[cfg(feature = "munsell")]
mod munsell_matt;

#[cfg(feature = "munsell")]
pub use munsell_matt::*;

use std::{
    borrow::Cow,
    ops::{Add, AddAssign, Deref, DerefMut, Mul, MulAssign},
};

use approx::{assert_abs_diff_eq, AbsDiffEq};
use nalgebra::SVector;

use crate::{
    error::Error,
    illuminant::CieIlluminant,
    lab::CieLab,
    physics::{gaussian_peak_one, wavelength},
    prelude::{Illuminant, Observer, D65, SPECTRUM_WAVELENGTH_RANGE},
    spectrum::{wavelengths, Spectrum, NS},
    traits::{Filter, Light},
};

/// # Colorant
///
/// A Colorant represents color filters and color patches, with spectral values between between 0.0
/// and 1.0, where 0.0 means all the light of a particular wavelength is absorbed, and 1.0 means no
/// light with that wavelength is absorbed, but reflected or transmitted instead.
///
/// Spectral values are represented as a vector of 401 values, covering the
/// wavelength range from 380 to 780 nanometer, including the end points, with a step size of 1
/// nanometer.  A spectral value in this vector has a unit of per nanometer. For example, a value of
/// 0.10 for a wavelength of 550 nanometer means that 10% of the light with that wavelength, a
/// greenish color, is not absorbed, but reflected or transmitted.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Colorant(pub(crate) Spectrum);

impl Colorant {
    /// Creates a Colorant from a spectrum, while validating the spectrum values.
    ///
    /// # Errors
    ///
    /// - CmtError::OutOfRange when the spectrum contains values outside the range 0.0 to 1.0.
    pub fn new(spectrum: Spectrum) -> Result<Self, Error> {
        if spectrum.values().iter().any(|v| !(0.0..=1.0).contains(v)) {
            Err(Error::OutOfRange {
                name: "Colorant Spectral Value".into(),
                low: 0.0,
                high: 1.0,
            })
        } else {
            Ok(Self(spectrum))
        }
    }

    /// Theoretical spectrum of a perfect grey colorant, consisting of 401
    /// values equal to the value given in the argument, over a range from 380
    /// to 780 nanometer. Mainly used for color mixing calculations.
    pub fn gray(gval: f64) -> Self {
        Self(Spectrum(SVector::<f64, NS>::repeat(gval.clamp(0.0, 1.0))))
    }

    /// Theoretical spectrum of a perfect white colorant, consisting of 401
    /// values over a range from 380 to 780 nanometer. Mainly used for
    /// color mixing calculations.
    pub fn white() -> Self {
        Self::gray(1.0)
    }

    /// Theoretical spectrum of a perfect black color patch, consisting of 401
    /// zero values over a range from 380 to 780 nanometer. Mainly used for
    /// color mixing calculations.
    pub fn black() -> Self {
        Self::gray(0.0)
    }

    /// A Rectangular Band Filter, specified by a central wavelength, and a
    /// width, both in units of meter, or nanometer.
    ///
    /// The filter has a peak value of 1.0
    /// ```rust
    /// # use approx::assert_ulps_eq;
    /// use colorimetry::prelude::*;
    /// let colorant = Colorant::top_hat(550.0, 1.0);
    /// let bandfilter = colorant.spectrum();
    /// assert_ulps_eq!(bandfilter[549], 0.0);
    /// assert_ulps_eq!(bandfilter[550], 1.0);
    /// assert_ulps_eq!(bandfilter[551], 0.0);
    ///
    /// let colorant = Colorant::top_hat(550.0, 2.0);
    /// let bandfilter = colorant.spectrum();
    /// assert_ulps_eq!(bandfilter[548], 0.0);
    /// assert_ulps_eq!(bandfilter[549], 1.0);
    /// assert_ulps_eq!(bandfilter[550], 1.0);
    /// assert_ulps_eq!(bandfilter[551], 1.0);
    /// assert_ulps_eq!(bandfilter[552], 0.0);
    ///
    /// ```
    pub fn top_hat(center: f64, width: f64) -> Self {
        let [center_m, width_m] = wavelengths([center, width]);
        let left = center_m - width_m / 2.0;
        let right = center_m + width_m / 2.0;
        let data = SVector::<f64, NS>::from_fn(|i, _j| {
            let w = wavelength(i + SPECTRUM_WAVELENGTH_RANGE.start());
            if w < left - f64::EPSILON || w > right + f64::EPSILON {
                0.0
            } else {
                1.0
            }
        });
        Self(Spectrum(data))
    }

    /// A Gaussian Filter, specified by a central wavelength, and a
    /// standard deviation `sigma` value, both in units of meter, or nanometer.
    ///
    /// The filter has a peak value of 1.0
    pub fn gaussian(center: f64, sigma: f64) -> Self {
        let [center_m, width_m] = wavelengths([center, sigma]);
        let data = SVector::<f64, NS>::from_fn(|i, _j| {
            gaussian_peak_one(
                wavelength(i + SPECTRUM_WAVELENGTH_RANGE.start()),
                center_m,
                width_m,
            )
        });
        Self(Spectrum(data))
    }

    /// Calculates the [Colorant] CIELAB values, using an illuminant and observer.
    ///
    /// The illuminant and observer are optional parameters.
    /// If illuminant is `None`, the D65 illuminant is used.
    /// If observer is `None`, the CIE 1931 observer is used.
    pub fn cielab(&self, illuminant_opt: Option<&dyn Light>, obs_opt: Option<Observer>) -> CieLab {
        let illuminant = illuminant_opt.unwrap_or(&D65);
        let obs = obs_opt.unwrap_or_default();
        let xyzn = obs.data().xyz(illuminant, None).set_illuminance(100.0);
        let xyz = obs.data().xyz(illuminant, Some(self));
        // unwrap is safe here, as we know the illuminant and observer are valid
        CieLab::from_xyz(xyz, xyzn).unwrap()
    }
}

#[test]
fn test_colorant_cielab() {
    // Test that the CIELAB values for a white colorant are as expected.
    // A white surface has CIELAB values of L* = 100, a* = 0, b* = 0.
    use crate::prelude::*;
    use approx::assert_abs_diff_eq;
    let colorant = Colorant::white();
    let [l, a, b] = colorant.cielab(None, None).values();
    assert_abs_diff_eq!(l, 100.0, epsilon = 1E-4); // L* should be 100 for white
    assert_abs_diff_eq!(a, 0.0, epsilon = 1E-4); // a* should be 0 for white
    assert_abs_diff_eq!(b, 0.0, epsilon = 1E-4); // b* should be 0 for white
}

impl TryFrom<Spectrum> for Colorant {
    type Error = Error;

    /// Creates a Colorant from a spectrum, while validating the spectrum values.
    ///
    /// # Errors
    ///
    /// Returns an error if the spectrum contains values outside the range 0.0 to 1.0.
    fn try_from(spectrum: Spectrum) -> Result<Self, Self::Error> {
        Self::new(spectrum)
    }
}

impl<F> From<F> for Colorant
where
    F: Fn(f64) -> f64,
{
    /**
        Colorant from an analytical function, defined over a domain from 0.0 to 1.0, covering the
        wavelength range from 380 to 780 nanometer.

        Values are clamped to a range from 0.0 to 1.0.
        ```rust
        use colorimetry::prelude::*;

        // linear filter from 0.0 to 1.0.
        let tilt: Colorant = (|x:f64|x).into();
        let xy = CIE1931.xyz(&D65, Some(&tilt)).chromaticity().to_array();
        approx::assert_abs_diff_eq!(xy.as_ref(), [0.4066, 0.4049].as_ref(), epsilon = 1E-4);

        // parabolic filter
        let parabolic: Colorant = (|x:f64|1.0 - 4.0 * (x - 0.5).powi(2)).into();
        let xy = CIE1931.xyz(&D65, Some(&parabolic)).chromaticity().to_array();
        approx::assert_abs_diff_eq!(xy.as_ref(), [0.3466, 0.3862].as_ref(), epsilon = 1E-4);
        ```
    */
    fn from(f: F) -> Self {
        let data = SVector::from_fn(|i, _j| {
            let x = i as f64 / (NS - 1) as f64;
            f(x).clamp(0.0, 1.0)
        });
        Colorant(Spectrum(data))
    }
}

/// Make colorant data available as a generic [`Filter`] entity, used in particular
/// in the [`Observer`] tristiumulus `xyz`-function.
impl Filter for Colorant {
    fn spectrum(&self) -> Cow<Spectrum> {
        Cow::Borrowed(&self.0)
    }
}

/// Adds together the spectrums of two colorants, resulting in a new colorant.
///
/// The result is clamped to the valid colorant range of 0.0 to 1.0.
impl Add for Colorant {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut spectrum = self.0 + rhs.0;
        spectrum.clamp(0.0, 1.0);
        Colorant(spectrum)
    }
}

/// Multiply the spectrum of a colorant with a scalar value.
///
/// The result is clamped to the valid colorant range of 0.0 to 1.0.
impl Mul<f64> for Colorant {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        let mut spectrum = self.0 * rhs;
        spectrum.clamp(0.0, 1.0);
        Self(spectrum)
    }
}

/// Multiply the spectrum of a colorant with a scalar value.
///
/// The result is clamped to the valid colorant range of 0.0 to 1.0.
impl Mul<Colorant> for f64 {
    type Output = Colorant;

    fn mul(self, rhs: Colorant) -> Self::Output {
        rhs.mul(self)
    }
}

impl Mul<Colorant> for Colorant {
    type Output = Self;

    /// Multiplication of two colorants using the `*`-operator.
    ///
    /// Subtractive Mixing.
    /// ```rust
    /// use colorimetry::prelude::*;
    /// let w = Colorant::white();
    /// let b = Colorant::black();
    /// let r: Colorant = w * b;
    /// ```
    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0) // use spectrum multiplication
    }
}

impl Mul<&Colorant> for &Colorant {
    type Output = Colorant;

    /// Multiplication of two colorant references using the `*`-operator.
    ///
    /// Non-consuming multiplication.
    /// Subtractive Mixing.
    /// ```rust
    /// use colorimetry::prelude::*;
    /// let w = Colorant::white();
    /// let b = Colorant::black();
    /// let r: Colorant = &w * &b;
    /// approx::assert_abs_diff_eq!(r,b);
    /// ```
    fn mul(self, rhs: &Colorant) -> Self::Output {
        Colorant(self.0 * rhs.0) // use spectrum multiplication
    }
}

/// Add the spectrum of a colorant to this colorant.
///
/// The result is clamped to the valid colorant range of 0.0 to 1.0.
impl AddAssign<&Self> for Colorant {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += rhs.0;
        self.0.clamp(0.0, 1.0);
    }
}

impl AbsDiffEq for Colorant {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.spectrum().abs_diff_eq(&other.spectrum(), epsilon)
    }
}
