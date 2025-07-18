//! # Spectral Power Distributions for Standard and User Definied Illuminants
//!
//! This module defines spectral **illuminants**—spectral power distributions (SPDs)
//! representing light sources or combinations of sources (e.g. daylight, LEDs, black-body radiators).
//! Spectral data cover 380 nm…780 nm at 1 nm intervals (401 samples).

mod cie_data;
pub use cie_data::*;

mod cie_illuminant;
pub use cie_illuminant::CieIlluminant;

#[cfg(feature = "cct")]
mod cct;

#[cfg(feature = "cct")]
pub use cct::CCT;

#[cfg(feature = "cri")]
mod cri;

mod led;
mod planck;
pub use planck::Planck;

#[cfg(feature = "cri")]
pub use cri::*;

#[cfg(feature = "cfi")]
mod cfi;

#[cfg(feature = "cfi")]
pub use cfi::CFI;

use std::{borrow::Cow, ops::Mul};

use nalgebra::{ArrayStorage, SMatrix, SVector};

use crate::{
    error::Error,
    observer::Observer,
    spectrum::{wavelength, wavelengths, Spectrum, NS, SPECTRUM_WAVELENGTH_RANGE},
    traits::Light,
    xyz::XYZ,
};

use led::led_ohno;

#[cfg(target_arch = "wasm32")]
mod wasm;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Clone, Default)]
/// # Illuminant
///
/// An illuminant is a spectral power distribution that represents the
/// spectral power density of a light source (sun, bulb, LED, etc.) in
/// W/m²/nm over 380–780 nm (401 samples).
pub struct Illuminant(pub(crate) Spectrum);

impl AsRef<Spectrum> for Illuminant {
    /// Allows using an `Illuminant` as a reference to its inner `Spectrum`.
    /// This is useful for passing an `Illuminant` to functions that expect a `&Spectrum`.
    ///
    /// # Examples
    /// ```rust
    /// use colorimetry::prelude::*;
    ///
    /// let illuminant = Illuminant::d65();
    /// let spectrum: &Spectrum = illuminant.as_ref();
    /// assert_eq!(spectrum.values().len(), 401);
    /// ```
    fn as_ref(&self) -> &Spectrum {
        &self.0
    }
}

impl Illuminant {
    /// Creates an illuminant directly from a spectrum.
    pub fn new(spectrum: Spectrum) -> Self {
        Illuminant(spectrum)
    }

    /// E, or Equal Energy Illuminant with an irradiance of 1 Watt per square
    /// meter in the spectrum between 380 and 780 nanometer
    pub fn equal_energy() -> Self {
        let s = 1. / NS as f64;
        Self(Spectrum(SVector::<f64, NS>::repeat(s)))
    }

    /**
    Standard Daylight Spectrum representing average daylight.

    It's truncated from the official standard, which
    covers 300 to 830 nanometers. It has a correlated color temperature of 6500
    K and should be used in color calculations requiring representative
    daylight. Variations occur based on factors like season, time of day, and
    location. For more details, refer to ISO 10526:1999/CIE
    S005/E-1998.
    */
    pub fn d65() -> Self {
        D65.clone()
    }

    /// CIE D50 Illuminant Standard Spectrum with 401 values over a range from
    /// 380 to 780 nanometers, with an interval size of 1 nanometer. Please be
    /// aware that this spectrum is truncated from the official standard, which
    /// is defined over a range from 300 to 830 nanometer.
    ///
    /// For most applications CIE recommends to use the D65 illuminant, to
    /// represent daylight, but this illuminant is often used in the printing
    /// industry.
    pub fn d50() -> Self {
        D50.clone()
    }

    /// A pure thermal emission based illuminant according to Planck's law.
    ///
    /// The generated spectrum is scaled to have a total power, over the full
    /// spectrum (including infrared), of 1 Watt.
    /// ```rust
    /// # use colorimetry::prelude::*;
    /// # use approx::assert_ulps_eq;
    ///
    /// let p3000 = Illuminant::planckian(3000.0);
    /// let xyz = Cie1931.xyz(&p3000, None);
    /// let chromaticity = xyz.chromaticity();
    /// assert_ulps_eq!(chromaticity.x(), 0.436_935, epsilon = 1E-6);
    /// assert_ulps_eq!(chromaticity.y(), 0.404_083, epsilon = 1E-6);
    ///
    /// ```
    pub fn planckian(cct: f64) -> Self {
        let p = Planck::new(cct);
        let s = 1E-9 / p.total_radiance(); // 1W/m2 total irradiance
        let data = SVector::<f64, NS>::from_fn(|i, _j| {
            s * p.at_wavelength((i + SPECTRUM_WAVELENGTH_RANGE.start()) as f64 * 1e-9)
        });
        Self(Spectrum(data))
    }

    /// Returns the reference illuminant for a test source with the given correlated color
    /// temperature (CCT) as defined by the TM-30-20 standard.
    ///
    /// The reference illuminant varies with the correlated color temperature (CCT)
    /// of the illuminant being tested. For a CCT below 4000 K, a blackbody Planckian
    /// illuminant is returned. For a CCT above 5000 K, a CIE D illuminant is returned.
    /// For CCTs between 4000 K and 5000 K, the two illuminants are blended to create a
    /// smooth crossover.
    pub fn cfi_reference(cct: f64) -> Result<Self, crate::Error> {
        const BLEND_RANGE_START: f64 = 4000.0;
        const BLEND_RANGE_END: f64 = 5000.0;

        if cct < BLEND_RANGE_START {
            Ok(Self::planckian(cct))
        } else if cct > BLEND_RANGE_END {
            Self::d_illuminant(cct)
        } else {
            let illuminant_planckian = Self::planckian(BLEND_RANGE_START);
            let illuminant_d = Self::d_illuminant(BLEND_RANGE_END)?;

            let ratio_d = (cct - BLEND_RANGE_START) / (BLEND_RANGE_END - BLEND_RANGE_START);
            let ratio_planckian = 1.0 - ratio_d;

            Ok(Self(
                ratio_d * illuminant_d.0 + ratio_planckian * illuminant_planckian.0,
            ))
        }
    }

    /// A spectral power distribution for a Light Emitting Diode.
    ///
    /// The spectrum is definded by a center wavelength, in units of meter or
    /// nanometer, and a full-width-half-maximum value, also in units of meter
    /// or nanometer. The generated spectrum is based on the model as published
    /// by Yoshi Ohno, from NIST, in his article, _Spectral Design
    /// considerations for white LED Color Rendering_, **Optical Engineering 44(11)**,
    /// November 2005.
    pub fn led(center: f64, width: f64) -> Self {
        let [center_m, width_m] = wavelengths([center, width]);
        let data = SVector::<f64, NS>::from_fn(|i,_j|
            // j = 0, first column
            led_ohno(wavelength(i + SPECTRUM_WAVELENGTH_RANGE.start()), center_m, width_m) * 1E-9);
        Self(Spectrum(data))
    }

    /// For an illuminant spectrum, the spectral values are scaled to have the specified
    /// irradiance, which is expressed in Watts per square meter.
    /// Typically, this is used to set the irradiance of an illuminant spectrum to 1.0 W/m².
    pub fn set_irradiance(mut self, irradiance: f64) -> Self {
        let s = irradiance / self.0 .0.sum();
        self.0 .0.iter_mut().for_each(|v| *v *= s);
        self
    }

    /// Calculate the irradiance of the illuminant spectrum, which is expressed in watts per square meter.
    pub fn irradiance(&self) -> f64 {
        self.0 .0.sum()
    }

    /// Sets the illuminance of the illuminant spectrum, which is expressed lumen per square meter,
    /// also referred to as lux.
    pub fn set_illuminance(mut self, obs: Observer, illuminance: f64) -> Self {
        let y = obs.y_from_spectrum(self.as_ref());
        let l = illuminance / y;
        self.0 .0.iter_mut().for_each(|v| *v *= l);
        self
    }

    /// Calculates the illuminance of the illuminant spectrum, which is expressed in lumen per square meter,
    /// also referred to as lux.
    pub fn illuminance(&self, obs: Observer) -> f64 {
        obs.y_from_spectrum(self.as_ref())
    }

    /// Calculates the Color Rendering Index values for illuminant spectrum.
    ///
    /// # Errors
    /// - CmtError::OutOfRange when the illuminant's distance to the Planckian locus is larger than 0.05 DUV,
    ///   or when the CCT is outside the range of 1000 to 25000 Kelvin.
    #[cfg(feature = "cri")]
    pub fn cri(&self) -> Result<CRI, Error> {
        self.try_into()
    }

    #[cfg(feature = "cfi")]
    pub fn cfi(&self) -> Result<CFI, Error> {
        CFI::new(self)
    }

    /// Creates a CIE D Illuminant with a correlated color temperature (CCT) in Kelvin.
    /// # Errors
    /// - CmtError::OutOfRange when the cct argument is below 4000 or above 25000 Kelvin.
    pub fn d_illuminant(cct: f64) -> Result<Illuminant, Error> {
        if !(4000.0..=25000.0).contains(&cct) {
            Err(Error::OutOfRange {
                name: "CIE D Illuminant Temperature".to_string(),
                low: 4000.0,
                high: 25000.0,
            })
        } else {
            let xd = match cct {
                t if t < 7000.0 => {
                    0.244063 + 0.09911E3 / t + 2.9678E6 / t.powi(2) - 4.607E9 / t.powi(3)
                }
                t => 0.23704 + 0.24748E3 / t + 1.9018E6 / t.powi(2) - 2.0064E9 / t.powi(3),
            };
            let yd = -3. * xd.powi(2) + 2.87 * xd - 0.275;
            let m = 0.0241 + 0.2562 * xd - 0.7341 * yd;
            let m1 = (-1.3515 - 1.7703 * xd + 5.9114 * yd) / m;
            let m2 = (0.03 - 31.4424 * xd + 30.0717 * yd) / m;
            let mut v = [0.0; CIE_D_S_LEN];
            v.iter_mut().enumerate().for_each(|(i, x)| {
                *x = CIE_D_S[(i, 0)] + m1 * CIE_D_S[(i, 1)] + m2 * CIE_D_S[(i, 2)]
            });
            let s = Spectrum::linear_interpolate(
                &[
                    *SPECTRUM_WAVELENGTH_RANGE.start() as f64,
                    *SPECTRUM_WAVELENGTH_RANGE.end() as f64,
                ],
                &v,
            )
            .unwrap();
            Ok(Illuminant(s).set_irradiance(1.0))
        }
    }

    /// Returns the XYZ tristimulus values for the illuminant.
    /// The values are calculated for the specified observer, or the default CIE 1931 observer if
    /// none is provided.
    pub fn xyz(&self, obs_opt: Option<Observer>) -> XYZ {
        let obs = obs_opt.unwrap_or_default();
        obs.xyz_from_spectrum(&self.0)
    }

    /// Calculate the correlated color temperature (CCT) of the illuminant.
    ///
    /// # Errors
    /// - CmtError::OutOfRange when the the distance to the Planckian locus is larger than 0.05 DUV,
    ///   or when the CCT is outside the range of 1000 to 25000 Kelvin.
    #[cfg(feature = "cct")]
    pub fn cct(&self) -> Result<cct::CCT, Error> {
        // CIE requires using the CIE1931 observer for calculating the CCT.
        let xyz = self.xyz(Some(Observer::Cie1931));
        cct::CCT::from_xyz(xyz)
    }
}

/// Creates an illuminant directly from a spectrum.
impl From<Spectrum> for Illuminant {
    fn from(spectrum: Spectrum) -> Self {
        Self::new(spectrum)
    }
}

impl Mul<f64> for Illuminant {
    /// Multiply a spectrum with a scalar f64 value.
    /// ```
    ///     use colorimetry::prelude::*;
    ///     use approx::assert_ulps_eq;
    ///
    ///     let mut led = Illuminant::led(550.0, 25.0);
    ///     let mut irradiance = led.irradiance();
    ///     assert_ulps_eq!(led.irradiance(), 1.0, epsilon = 1E-10);
    ///
    ///     led = led * 10.0;
    ///     assert_ulps_eq!(led.irradiance(), 10.0, epsilon = 1E-10);
    /// ```
    type Output = Self;

    // spectrum * scalar
    fn mul(self, rhs: f64) -> Self::Output {
        Self(self.0 * rhs) // uses Spectrum::Mul
    }
}

impl Mul<Illuminant> for f64 {
    /// Multiply a spectrum with a scalar f64 value.
    /// ```
    ///     use colorimetry::prelude::*;
    ///     use approx::assert_ulps_eq;
    ///
    ///     let mut led = Illuminant::led(550.0, 25.0);
    ///     let mut irradiance = led.irradiance();
    ///     assert_ulps_eq!(led.irradiance(), 1.0, epsilon = 1E-10);
    ///
    ///     led = 10.0 * led;
    ///     assert_ulps_eq!(led.irradiance(), 10.0, epsilon = 1E-10);
    /// ```
    type Output = Illuminant;

    // scalar * spectrum
    fn mul(self, rhs: Illuminant) -> Self::Output {
        Illuminant(self * rhs.0)
    }
}

impl Light for Illuminant {
    fn spectrum(&self) -> Cow<'_, Spectrum> {
        Cow::Borrowed(self.as_ref())
    }
}

#[test]
fn test_d_illuminant() {
    use crate::observer::Observer::Cie1931;
    use crate::prelude::*;
    let s = Illuminant::d_illuminant(6504.0).unwrap();
    let xyz = Cie1931.xyz_from_spectrum(s.as_ref()).set_illuminance(100.0);
    approx::assert_ulps_eq!(xyz, Cie1931.xyz_d65(), epsilon = 2E-2);
}

#[test]
fn test_d_illuminant_range_error() {
    use crate::prelude::*;
    let s = Illuminant::d_illuminant(3999.0);
    assert!(s.is_err());
    let s = Illuminant::d_illuminant(25001.0);
    assert!(s.is_err());
}

#[test]
fn test_xyz() {
    use crate::prelude::*;
    let s = *Illuminant::d_illuminant(6504.0).unwrap().as_ref().values();
    let illuminant = Illuminant(Spectrum::from(s));
    let xyz = illuminant.xyz(None).set_illuminance(100.0);
    approx::assert_ulps_eq!(xyz, Observer::Cie1931.xyz_d65(), epsilon = 2E-2);
}

const CIE_D_S_LEN: usize = 81;

static CIE_D_S: SMatrix<f64, CIE_D_S_LEN, 3> = SMatrix::from_array_storage(ArrayStorage([
    [
        63.40, 64.60, 65.80, 80.30, 94.80, 99.80, 104.80, 105.35, 105.90, 101.35, 96.80, 105.35,
        113.90, 119.75, 125.60, 125.55, 125.50, 123.40, 121.30, 121.30, 121.30, 117.40, 113.50,
        113.30, 113.10, 111.95, 110.80, 108.65, 106.50, 107.65, 108.80, 107.05, 105.30, 104.85,
        104.40, 102.20, 100.00, 98.00, 96.00, 95.55, 95.10, 92.10, 89.10, 89.80, 90.50, 90.40,
        90.30, 89.35, 88.40, 86.20, 84.00, 84.55, 85.10, 83.50, 81.90, 82.25, 82.60, 83.75, 84.90,
        83.10, 81.30, 76.60, 71.90, 73.10, 74.30, 75.35, 76.40, 69.85, 63.30, 67.50, 71.70, 74.35,
        77.00, 71.10, 65.20, 56.45, 47.70, 58.15, 68.60, 66.80, 65.00,
    ],
    [
        38.50, 36.75, 35.00, 39.20, 43.40, 44.85, 46.30, 45.10, 43.90, 40.50, 37.10, 36.90, 36.70,
        36.30, 35.90, 34.25, 32.60, 30.25, 27.90, 26.10, 24.30, 22.20, 20.10, 18.15, 16.20, 14.70,
        13.20, 10.90, 8.60, 7.35, 6.10, 5.15, 4.20, 3.05, 1.90, 0.95, 0.00, -0.80, -1.60, -2.55,
        -3.50, -3.50, -3.50, -4.65, -5.80, -6.50, -7.20, -7.90, -8.60, -9.05, -9.50, -10.20,
        -10.90, -10.80, -10.70, -11.35, -12.00, -13.00, -14.00, -13.80, -13.60, -12.80, -12.00,
        -12.65, -13.30, -13.10, -12.90, -11.75, -10.60, -11.10, -11.60, -11.90, -12.20, -11.20,
        -10.20, -9.00, -7.80, -9.50, -11.20, -10.80, -10.40,
    ],
    [
        3.00, 2.10, 1.20, 0.05, -1.10, -0.80, -0.50, -0.60, -0.70, -0.95, -1.20, -1.90, -2.60,
        -2.75, -2.90, -2.85, -2.80, -2.70, -2.60, -2.60, -2.60, -2.20, -1.80, -1.65, -1.50, -1.40,
        -1.30, -1.25, -1.20, -1.10, -1.00, -0.75, -0.50, -0.40, -0.30, -0.15, 0.00, 0.10, 0.20,
        0.35, 0.50, 1.30, 2.10, 2.65, 3.20, 3.65, 4.10, 4.40, 4.70, 4.90, 5.10, 5.90, 6.70, 7.00,
        7.30, 7.95, 8.60, 9.20, 9.80, 10.00, 10.20, 9.25, 8.30, 8.95, 9.60, 9.05, 8.50, 7.75, 7.00,
        7.30, 7.60, 7.80, 8.00, 7.35, 6.70, 5.95, 5.20, 6.30, 7.40, 7.10, 6.80,
    ],
]));
