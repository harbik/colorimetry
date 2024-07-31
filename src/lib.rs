
//! # Overview
//!
//!
//! `colorimetry` is a Rust and WebAssembly/JavaScript toolbox for color calculations in illumination and engineering projects.
//! It implements algorithms as recommended by the International Commission on Illumination,
//! the International Color Consortium, the Illumination Engineering Society, and many others.
//!
//! 
//! This example calculates the Illuminance and CIE 1931 (x, y) chromaticity
//! coordinates for a Planckian (thermal emission based) illuminator, with a
//! Correlated Color Temperature of 3000 Kelvin.
//!
//! ```rust
//! use crate::colorimetry::{Spectrum, CIE1931};
//! use approx::assert_ulps_eq;
//! 
//! let p3000 = Spectrum::planckian(3000.0);
//! let [l, x, y] = CIE1931.xyz(&p3000).lxy();
//!
//! assert_ulps_eq!(l, 20.668_927, epsilon = 1E-6);
//! assert_ulps_eq!(x, 0.436_935, epsilon = 1E-6);
//! assert_ulps_eq!(y, 0.404_083, epsilon = 1E-6);
//! ```
//! 
//! The `Spectrum` class has many constructors, and in this case the `planckian` constructor is used.
//! Many other illuminant constructors are availables, such as CIE D65, D50, A, Led.
//! 
//! `Spectrum` data are defined over a domain from 380 to 780 nanometer, with an interval size of 1 nanometer.
//! For Illuminants the values are supposed to have an SI unit of Watt per square meter per nanometer.
//! 
//! `CIE1931` is and instance of the `Observer` class representing colorimetric
//! standard observers, and also includes the CIE 1976 10º standard observers,
//! and the CIE 2015 2º and 10º cone fundamental derived observers.
//! It's main function `xyz(s: &Spectrum) -> XYZ`, which takes a spectral distribution as single argument,
//! with an instance of the `XYZ` class, which encapsulates the X, Y, and Z tristimulus values.
//!
//! Tristimulus values are the basis for many color madels.

#![allow(dead_code)]

pub use spc::Spectrum;
pub use obs::ObsId;
pub use data::CIE1931;


pub mod obs;
pub mod spc;
pub mod xyz;
pub mod lab;
pub mod rgb;
pub mod physics;
pub mod data;