use std::ops::Deref;

use colored::Color;

use crate::Spectrum;


pub struct Colorant(Spectrum);

impl Deref for Colorant {
    type Target = Spectrum;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}