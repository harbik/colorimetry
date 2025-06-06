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
//!   - Standard observers CIE1931
//! - **Physics**  
//!   - Radiometric constants (Planck’s law, Stefan–Boltzmann), Gaussian helpers, wavelength converters.  
//! - **RGB Color**  
//!   - `Rgb`, `WideRgb`, color spaces, gamma curves, and gamut-clamped conversions.  
//! - **Spectra & Stimuli**  
//!   - `Spectrum`, `Stimulus`, wavelength-domain utilities, and prelude traits (`Light`, `Filter`).  
//! - **XYZ Tristimulus**  
//!   - The `XYZ` type, chromaticity, and conversion helpers.  
//!

pub use super::cam::{CamTransforms, CieCam16, ViewConditions};
pub use super::colorant::Colorant;
pub use super::illuminant::{CieIlluminant, Illuminant};
pub use super::lab::CieLab;
pub use super::observer::{Observer, Observer::Cie1931};
pub use super::rgb::{Rgb, RgbSpace, WideRgb};
pub use super::spectrum::Spectrum;
pub use super::stimulus::Stimulus;
pub use super::traits::{Filter, Light};
pub use super::xyz::{Chromaticity, XYZ};
