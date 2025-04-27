use std::{
    borrow::Cow,
    ops::{Add, AddAssign, Deref, DerefMut, Mul, MulAssign},
};

use approx::{assert_abs_diff_eq, AbsDiffEq};
use nalgebra::SVector;

use crate::{
    error::CmtError,
    lab::CieLab,
    physics::{gaussian_peak_one, wavelength},
    prelude::{Illuminant, Observer, D65},
    spectrum::{wavelengths, Spectrum, NS},
    std_illuminants,
    traits::{Filter, Light},
};

/// A Colorant represents color filters and color patches. These are spectrums
/// with spectral values between 0.0 and 1.0.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Colorant(pub(crate) Spectrum);

impl Colorant {
    /// Creates a Colorant from a spectrum, while validating the spectrum values.
    ///
    /// # Errors
    ///
    /// Returns an error if the spectrum contains values outside the range 0.0 to 1.0.
    pub fn new(spectrum: Spectrum) -> Result<Self, CmtError> {
        if spectrum.values().iter().any(|v| !(0.0..=1.0).contains(v)) {
            Err(CmtError::OutOfRange {
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
            //  let w = (i+380) as f64 * 1E-9;
            let w = wavelength(i + 380);
            if w < left - f64::EPSILON || w > right + f64::EPSILON {
                0.0
            } else {
                1.0
            }
        });
        Self(Spectrum(data))
    }

    /// A Gaussian Filter, specified by a central wavelength, and a
    /// full-width-half-maximum value, both in units of meter, or nanometer.
    ///
    /// The filter has a peak value of 1.0
    pub fn gaussian(center: f64, sigma: f64) -> Self {
        let [center_m, width_m] = wavelengths([center, sigma]);
        let data = SVector::<f64, NS>::from_fn(|i, _j| {
            gaussian_peak_one((i + 380) as f64 * 1E-9, center_m, width_m)
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
        let xyz = obs
            .data()
            .xyz(illuminant, Some(self))
            .set_illuminance(100.0);
        CieLab::try_from(xyz).unwrap()
    }
}

#[test]
fn test_colorant_cielab() {
    // Test that the CIELAB values for a white colorant are as expected.
    // A white surface has CIELAB values of L* = 100, a* = 0, b* = 0.
    use crate::prelude::*;
    use approx::assert_abs_diff_eq;
    let colorant = Colorant::white();
    let lab = colorant.cielab(None, None);
    assert_abs_diff_eq!(lab.lab[0], 100.0, epsilon = 1E-4); // L* should be 100 for white
    assert_abs_diff_eq!(lab.lab[1], 0.0, epsilon = 1E-4); // a* should be 0 for white
    assert_abs_diff_eq!(lab.lab[2], 0.0, epsilon = 1E-4); // b* should be 0 for white
}

impl TryFrom<Spectrum> for Colorant {
    type Error = CmtError;

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
        let xy = CIE1931.xyz(&D65, Some(&tilt)).chromaticity();
        approx::assert_abs_diff_eq!(xy.as_ref(), [0.4066, 0.4049].as_ref(), epsilon = 1E-4);

        // parabolic filter
        let parabolic: Colorant = (|x:f64|1.0 - 4.0 * (x - 0.5).powi(2)).into();
        let xy = CIE1931.xyz(&D65, Some(&parabolic)).chromaticity();
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
/// in the [`Observer`](crate::observer::Observer) tristiumulus `xyz`-function.
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
