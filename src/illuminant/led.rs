// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

//! # Functions and Constants from Physics
//!
//! This module provides core physical constants and spectral‐generation functions
//! used throughout the library to build illuminants, LEDs, black‐body radiators,
//! and band-pass filters.  
//!
//! ## Spectral Band & LED Models
//! - `sigma_from_fwhm(fwhm)` / `fwhm_from_sigma(σ)` — Convert between FWHM and standard deviation.  
//! - `gaussian_peak_one(x, μ, σ)` — Un‐normalized Gaussian with peak = 1.0.  
//! - `gaussian_normalized(x, μ, σ)` — Standard Gaussian PDF.  
//! - `led_ohno(λ, center, width)` — NIST Ohno white‐LED SPD model (Optical Engineering 44(11), 2005).  
//!

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
