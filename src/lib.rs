
#![allow(dead_code)]
#![doc = include_str!("../README.md")]


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