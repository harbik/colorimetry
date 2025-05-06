#![allow(dead_code, unused_variables, unused_imports)]
#![doc = include_str!("../README.md")]
// This library defines many floating-point constants for colorimetry.
// Clippy's `approx_constant` lint would otherwise generate numerous false positives
// by flagging constants close to standard values.
// To suppress these, `#![allow(clippy::approx_constant)]` is applied.
#![allow(clippy::approx_constant)]

pub mod cam;
#[cfg(feature = "cct")]
pub mod cct;
pub mod colorant;
#[cfg(feature = "cri")]
pub mod cri;
pub mod data;
pub mod error;
pub mod gamma;
pub mod geometry;
pub mod illuminant;
pub mod lab;
#[cfg(feature = "munsell")]
pub mod munsell_matt;
pub mod observer;
pub mod physics;
pub mod prelude;
pub mod rgb;
pub mod rgbspace;
pub mod spectrum;
pub mod std_illuminants;
pub mod stimulus;
pub mod traits;
pub mod viewconditions;
pub mod xyz;
