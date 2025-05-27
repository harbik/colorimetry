#![allow(dead_code, unused_variables, unused_imports)]
#![doc = include_str!("../README.md")]
// This library defines many floating-point constants for colorimetry.
// Clippy's `approx_constant` lint would otherwise generate numerous false positives
// by flagging constants close to standard values.
// To suppress these, `#![allow(clippy::approx_constant)]` is applied.
#![allow(clippy::approx_constant)]

pub mod cam;
pub mod colorant;
pub mod error;
pub mod geometry;
pub mod illuminant;
pub mod lab;
pub mod observer;
pub mod physics;
pub mod rgb;
pub mod spectrum;
pub mod stimulus;
pub mod traits;
pub mod xyz;
