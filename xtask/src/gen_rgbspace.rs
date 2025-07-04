#![allow(unused)]
use std::{collections::HashMap, hash::Hash};

use colorimetry::{observer::Observer, prelude::CieIlluminant, rgb::GammaCurve, spectrum::NS, xyz::XYZ};
use strum::IntoEnumIterator;


pub struct GenRgbSpaceData {
    name: String,

    red_chromaticity_cie1931: [f64; 2],
    green_chromaticity_cie1931: [f64; 2],
    blue_chromaticity_cie1931: [f64; 2],

    red_stimulus: [f64; NS],
    green_stimulus: [f64; NS],
    blue_stimulus: [f64; NS],
    
    white: CieIlluminant,
    white_points: HashMap<Observer, XYZ>,
    
    gamma: GammaCurve,
}

impl GenRgbSpaceData {
    pub fn srgb() -> Self {
        todo!()
    }

    pub fn adobe_rgb() -> Self {
        todo!()
    }

    pub fn display_p3() -> Self {
        todo!()
    }

    pub fn cie_rgb() -> Self {
        todo!()
    }

    pub fn dci_p3() -> Self {
        todo!()
    }

    pub fn rec709() -> Self {
        todo!()
    }

    pub fn rec2020() -> Self {
        todo!()
    }

    pub fn prophoto_rgb() -> Self {
        todo!()
    }

    pub fn aces() -> Self {
        todo!()
    }

    pub fn acescg() -> Self {
        todo!()
    }

}



pub fn main() {
    // for each observer
    for observer in Observer::iter() {
        println!("{}", observer);
    }
}
