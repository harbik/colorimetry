#![allow(dead_code, unused_variables, unused_imports)]
#![doc = include_str!("../README.md")]
// This library contains lots and lots of float constants for colorimetry.
// If clippy finds a constant that is close to one of the well known standard library
// constanst, it will complain by default. This yields a lot of false positives.
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

// Set "rust-analyzer.check.features": "all" or ["cri", ...] to limit processing time
// `cargo build --all-features`
// `wasm-pack build --target web --release --all-features`
// or `cargo build --no-default-features`: 193K wasm-file
// `wasm-pack build --target web --release --no-default-features`
