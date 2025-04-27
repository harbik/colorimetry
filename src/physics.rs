use num_traits::ToPrimitive;
use wasm_bindgen::prelude::wasm_bindgen;

/// The speed of light (m/s)
pub const C: f64 = 299792458.0;

/// Boltzmann constant (m<sup>2</sup> kg s<sup>-2</sup> K<sup>-1</sup>)
pub const KB: f64 = 1.3806485279E-23;

/// Planck constant (m<sup>2</sup> kg / s)
pub const H: f64 = 6.6260700408181E-34;

/// First radiation constant (W m<sup>2</sup>)
pub const C1: f64 = 2. * std::f64::consts::PI * H * C * C;

/// Second radiation constant \( c_2 \) appears in Planck's law
/// and has the SI unit **m·K** (meter times kelvin)
/// It's definition is now exact, and is used to define the temperature scale, but it's value varied in the past.
/// see [Wiki ITS](https://en.wikipedia.org/wiki/Planckian_locus#International_Temperature_Scale)
pub const C2: f64 = H * C / KB;

/// Value as used in the definition oof the A Illuminant
pub const C2_NBS_1931: f64 = 1.435E-2;

/// Value as used in the D Illuminant series.
pub const C2_IPTS_1948: f64 = 1.4380E-2;
pub const C2_ITS_1968: f64 = 1.4388E-2;
//pub const FWHM2STDDEV: f64 = 2.354820045;

// calculated on first dereference, can not use floating point calculations in const (yet?)
pub(crate) static FWHM: LazyLock<f64> = LazyLock::new(|| (8.0 * 2f64.ln()).sqrt());

/**
Planck with the second radiant constant as parameter.

This to be used for planckian radiators not in vacuum, or to calculate standard lights defined with
older values of this constant.

l: wavelength in meter
t: absolute temperature in Kelvin
c2: second radiative constant, in meter * Kelvin; can also be used to include refractive index, using c2 ::  c2 / n

*/
#[inline]
pub fn planck_c2(l: f64, t: f64, c2: f64) -> f64 {
    crate::physics::C1 / l.powi(5) / ((c2 / (l * t)).exp() - 1.0)
}

/// Planck Temperature derivate: d(Planck)/dT
pub fn planck_slope_c2(l: f64, t: f64, c2: f64) -> f64 {
    let c3 = C1 * c2 / t.powi(2);
    let e = (c2 / (l * t)).exp();
    c3 / l.powi(6) * e / (e - 1.0).powi(2)
}

/// Planck Temperature second derivative: d2(Planck)/dT2
pub fn planck_curvature_c2(l: f64, t: f64, c2: f64) -> f64 {
    let e = (c2 / (l * t)).exp();
    planck_slope_c2(l, t, c2) / t * (c2 / (l * t) * (e + 1.0) / (e - 1.0) - 2.0)
}

#[inline]
pub fn planck(l: f64, t: f64) -> f64 {
    planck_c2(l, t, crate::physics::C2)
}

#[inline]
pub fn planck_slope(l: f64, t: f64) -> f64 {
    planck_slope_c2(l, t, crate::physics::C2)
}

/// Stefan-Boltzmann constant (W m<sup>-2</sup> K<sup>-4</sup>)
const SIGMA: f64 = 5.670_374_419_184E-8;

/// Stefan Boltzmann law: Blackbody's radiant emittance (W m<sup>-2</sup>), as function of its absolute
/// temperature (K).
#[inline]
#[wasm_bindgen(js_name= stefanBoltzmann)]
pub fn stefan_boltzmann(temperature: f64) -> f64 {
    SIGMA * temperature.powi(4)
}

const A: f64 = 1.11926158998; // scaling factor for power
const B: f64 = 1.08480681239; // scaling factor for width (l_w = B * l_fwhm)

/// LED Spectrum model
///
/// See Ohno, Spectral Design considerations for white LED Color Rendering, Optical Engineering 44(11), November 2005
/// Scale by spectralWidth
pub fn led_ohno(wl: f64, center: f64, width: f64) -> f64 {
    let width = B * width;
    let t = (wl - center) / (width);
    let g = libm::expm1(-(t.powi(2))) + 1.0;
    (g + 2.0 * g.powi(5)) / (3.0 * A * width)
}

use core::f64;
use std::{f64::consts::PI, sync::LazyLock};

#[inline]
pub fn sigma_from_fwhm(fwhm: f64) -> f64 {
    fwhm / *FWHM
}

#[inline]
pub fn fwhm_from_sigma(sigma: f64) -> f64 {
    sigma * *FWHM
}

#[inline]
pub fn gaussian_peak_one(x: f64, mu: f64, sigma: f64) -> f64 {
    let exponent = -((x - mu).powi(2)) / (2.0 * sigma.powi(2));
    exponent.exp()
}

#[test]
fn gaussian_peak_one_test() {
    use approx::assert_ulps_eq;
    let sigma = 10E-9;
    let mu = 500E-9;
    let x = mu - sigma;
    let v = gaussian_peak_one(x, mu, sigma);
    assert_ulps_eq!(v, 0.60653065971, epsilon = 1E-10);
}

#[inline]
pub fn gaussian_normalized(x: f64, mu: f64, sigma: f64) -> f64 {
    let exponent = -((x - mu).powi(2)) / (2.0 * sigma.powi(2));
    (1.0 / (sigma * (2.0 * PI).sqrt())) * exponent.exp()
}

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
/// use colorimetry::prelude::*;
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
    380E-9 * (1.0 - f) + 780E-9 * f
}
