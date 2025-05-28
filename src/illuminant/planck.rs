//! # Functions and Constants from Physics
//!
//! This module provides core physical constants and spectral‐generation functions
//! used throughout the library to build illuminants, LEDs, black‐body radiators,
//! and band-pass filters.  
//!
//! ## Constants
//! - `C` — Speed of light (m/s)  
//! - `KB` — Boltzmann constant (m²·kg·s⁻²·K⁻¹)  
//! - `H` — Planck constant (m²·kg/s)  
//! - `C1` — First radiation constant (W·m²)  
//! - `C2`, `C2_NBS_1931`, `C2_IPTS_1948`, `C2_ITS_1968` — Second radiation constants for various standards  
//! - `SIGMA` — Stefan–Boltzmann constant (W·m⁻²·K⁻⁴)  
//!
//! ## Black‐Body Radiation
//! - `planck_c2(l, T, c2)` — Planck’s law at wavelength `l` (m), temperature `T` (K), using constant `c2`.  
//! - `planck_slope_c2`, `planck_curvature_c2` — First and second derivatives with respect to temperature.  
//! - `planck(l, T)` / `planck_slope(l, T)` — Convenience wrappers using the current `C2`.  
//! - `stefan_boltzmann(T)` — Total radiant emittance of a black‐body at `T` (W·m⁻²).  
//!


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
    C1 / l.powi(5) / ((c2 / (l * t)).exp() - 1.0)
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
    planck_c2(l, t, C2)
}

#[inline]
pub fn planck_slope(l: f64, t: f64) -> f64 {
    planck_slope_c2(l, t, C2)
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
