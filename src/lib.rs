#![allow(dead_code, unused_variables)]
#![doc = include_str!("../README.md")]
// This library defines many floating-point constants for colorimetry.
// Clippy's `approx_constant` lint would otherwise generate numerous false positives
// by flagging constants close to standard values.
// To suppress these, `#![allow(clippy::approx_constant)]` is applied.
#![allow(clippy::approx_constant)]
///#![allow(dead_code, unused_variables, unused_imports)]
pub mod cam;
pub mod colorant;
mod error;
pub mod geometry;
pub mod illuminant;
pub mod lab;
pub mod math;
pub mod observer;
pub mod prelude;
pub mod rgb;
pub mod spectrum;
pub mod stimulus;
pub mod traits;
pub mod xyz;

pub use error::Error;
