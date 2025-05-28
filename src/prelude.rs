//! # Single Import for Colorimetry Types and Functions
//!
//! The `prelude` module gathers and re-exports the most commonly used types, traits, and functions
//! from across the Colorimetry library, so you can bring everything you need into scope with a
//! single import:
//!
//! ```rust
//! use colorimetry::prelude::*;
//! ```
//!
//! ## What’s in the Prelude
//!
//! - **Color Appearance**  
//!   - `CieCam16`, `ViewConditions`, and related CAM utilities.  
//! - **Colorants & Munsell**  
//!   - `Colorant`, plus `munsell_matt` types when the “munsell” feature is enabled.  
//! - **Geometry**  
//!   - Basic geometric utilities used in color calculations.  
//! - **Illuminants**  
//!   - Standard illuminants (`D65`, `D50`, etc.), Planckian and LED models, plus CRI/CCT support if enabled.  
//! - **CIELab & Delta-E**  
//!   - `CieLab` and color-difference formulas (Euclidean ΔE*ab, CIEDE2000).  
//! - **Observers**  
//!   - Standard observers (CIE1931, CIE1964, CIE2015 cone fundamentals).  
//! - **Physics**  
//!   - Radiometric constants (Planck’s law, Stefan–Boltzmann), Gaussian helpers, wavelength converters.  
//! - **RGB Color**  
//!   - `Rgb`, `WideRgb`, color spaces, gamma curves, and gamut-clamped conversions.  
//! - **Spectra & Stimuli**  
//!   - `Spectrum`, `Stimulus`, wavelength-domain utilities, and prelude traits (`Light`, `Filter`).  
//! - **XYZ Tristimulus**  
//!   - The `XYZ` type, chromaticity, and conversion helpers.  
//!
//! By importing `prelude::*`, you get a comprehensive set of tools for colorimetry tasks—no more
//! hunting through modules for the right name!

pub use super::cam::*;
pub use super::colorant::*;
pub use super::geometry::*;
pub use super::illuminant::*;
pub use super::lab::*;
pub use super::observer::*;
pub use super::physics::*;
pub use super::rgb::*;
pub use super::spectrum::*;
pub use super::stimulus::*;
pub use super::traits::*;
pub use super::xyz::*;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::JsValue;
