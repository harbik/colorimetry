//! Provides functionality for handling relative XYZ color space gamuts with different observers and illuminants.
//!
//! This module implements gamut handling for various CIE standard observers and illuminants, allowing for:
//! - Chromaticity binning and conversion
//! - Maximum luminance determination for specific chromaticity coordinates
//! - Conversion between bin coordinates and RelXYZ values
//!
//! # Features
//! - Default support for CIE 1931 2° observer with D50 and D65 illuminants
//! - Optional support for additional observers (CIE 1964, CIE 2015) via "supplemental-observers" feature
//! - High-resolution chromaticity binning (2000x2000 bins)

use std::collections::HashMap;

use crate::illuminant::CieIlluminant;
use crate::lab::CieLab;
use crate::observer::{Observer, OptimalColors};
use crate::traits::Light;
use crate::xyz::{RelXYZ, XYZ};

pub struct CieLChGamut {
    illuminant: CieIlluminant,
    observer: Observer,
    max_chromas: HashMap<[u16; 2], u16>,
    white_point: XYZ,
}

impl CieLChGamut {
    const H_BINS: u16 = 72;
    const L_BINS: u16 = 100;
    const H_SCALE: f64 = Self::H_BINS as f64 / 360.0;
    const L_SCALE: f64 = Self::L_BINS as f64 / 100.0;

    pub fn new(observer: Observer, illuminant: CieIlluminant) -> Self {
        let opt_colors = OptimalColors::new(observer, illuminant);
        let whitepoint = illuminant.xyzn100(observer);
        let mut max_chromas = HashMap::new();
        for xyz in opt_colors.colors().iter() {
            let rel_xyz = RelXYZ::from_vec(*xyz, opt_colors.white_point());
            let [l, c, h] = CieLab::from_xyz(rel_xyz).lch();

            let l_bin = Self::l_to_bin(l);
            let h_bin = Self::h_to_bin(h);
            let c_bin = Self::c_to_bin(c);

            // Only add if the discrete ab coordinates and luminance are valid.
            // This does XYZ - CieLab round trip, which is not ideal,
            // but it is a reasonable approximation for the purpose of this mapping.
            let cielab_for_bin = Self::bins_to_cielab_static(whitepoint, l_bin, c_bin, h_bin);
            if cielab_for_bin.is_valid() {
                // Insert or update the maximum luminance for this chromaticity bin
                max_chromas
                    .entry([l_bin, h_bin])
                    .and_modify(|existing| {
                        if *existing < c_bin {
                            *existing = c_bin;
                        }
                    })
                    .or_insert(c_bin);
            }
        }
        CieLChGamut {
            illuminant,
            observer,
            max_chromas,
            white_point: illuminant.xyzn100(observer),
        }
    }

    pub fn illuminant(&self) -> CieIlluminant {
        self.illuminant
    }

    pub fn observer(&self) -> Observer {
        self.observer
    }

    pub fn white_point(&self) -> XYZ {
        self.white_point
    }

    pub fn bins_to_cielab_static(whitepoint: XYZ, l_bin: u16, c_bin: u16, h_bin: u16) -> CieLab {
        let l = Self::bin_to_l(l_bin);
        let c = Self::bin_to_c(c_bin);
        let h = Self::bin_to_h(h_bin);
        CieLab::from_lch([l, c, h], whitepoint)
    }

    pub fn h_to_bin(value: f64) -> u16 {
        (value * Self::H_SCALE).floor() as u16
    }

    pub fn bin_to_h(bin: u16) -> f64 {
        (bin.min(Self::H_BINS) as f64) / Self::H_BINS as f64 * 360.0
    }

    pub fn c_to_bin(value: f64) -> u16 {
        value.floor() as u16
    }

    pub fn bin_to_c(bin: u16) -> f64 {
        bin as f64
    }

    pub fn l_to_bin(value: f64) -> u16 {
        (value * Self::L_SCALE).floor() as u16
    }

    pub fn bin_to_l(bin: u16) -> f64 {
        (bin.min(Self::L_BINS) as f64) / Self::L_BINS as f64 * 100.0
    }

    pub fn bins_to_cielab(&self, l_bin: u16, c_bin: u16, h_bin: u16) -> CieLab {
        Self::bins_to_cielab_static(self.white_point, l_bin, c_bin, h_bin)
    }

    pub fn max_chroma(&self, x: u16, y: u16) -> Option<u16> {
        self.max_chromas.get(&[x, y]).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_relxyz_gamut_data() {
        let gamut = CieLChGamut::new(Observer::Cie1931, CieIlluminant::D65);
        assert_eq!(gamut.observer(), Observer::Cie1931);
        assert!(!gamut.max_chromas.is_empty());
        let size_guess = (CieLChGamut::H_BINS * CieLChGamut::L_BINS) / 2;
        assert!(
            gamut.max_chromas.len() > size_guess as usize,
            "Size should be at least {}, got {}",
            size_guess,
            gamut.max_chromas.len()
        );
    }

    #[test]
    fn test_max_chroma_for_bin() {
        let gamut = CieLChGamut::new(Observer::Cie1931, CieIlluminant::D65);
        let l = 50;
        for h in 0..=100 {
            if let Some(c) = gamut.max_chroma(l, h) {
                println!("Max luminance found for bin ({l}, {h}): {c}");
                return;
            }
        }
        panic!("No max luminance found for any bin in the range 200-300");
    }

    #[test]
    fn test_bin_conversions() {
        // Test L* conversions
        assert_eq!(CieLChGamut::l_to_bin(50.0), 50);
        assert_abs_diff_eq!(CieLChGamut::bin_to_l(50), 50.0, epsilon = 0.001);
        assert_eq!(CieLChGamut::l_to_bin(100.0), 100);
        assert_abs_diff_eq!(CieLChGamut::bin_to_l(100), 100.0, epsilon = 0.001);

        // Test h conversions
        assert_eq!(CieLChGamut::h_to_bin(180.0), 36);
        assert_abs_diff_eq!(CieLChGamut::bin_to_h(36), 180.0, epsilon = 0.001);
        assert_eq!(CieLChGamut::h_to_bin(360.0), 72);
        assert_abs_diff_eq!(CieLChGamut::bin_to_h(72), 360.0, epsilon = 0.001);

        // Test C* conversions
        assert_eq!(CieLChGamut::c_to_bin(50.5), 50);
        assert_abs_diff_eq!(CieLChGamut::bin_to_c(50), 50.0, epsilon = 0.001);
    }

    #[test]
    fn test_white_point() {
        let gamut = CieLChGamut::new(Observer::Cie1931, CieIlluminant::D65);
        let wp = gamut.white_point();
        assert_abs_diff_eq!(wp.x(), 95.0422, epsilon = 0.001);
        assert_abs_diff_eq!(wp.y(), 100.0, epsilon = 0.001);
        assert_abs_diff_eq!(wp.z(), 108.8610, epsilon = 0.001);
    }

    #[test]
    fn test_bins_to_cielab() {
        let gamut = CieLChGamut::new(Observer::Cie1931, CieIlluminant::D65);
        let lab = gamut.bins_to_cielab(50, 10, 36);
        assert!(lab.is_valid());
        let [l, c, h] = lab.lch();
        assert_abs_diff_eq!(l, 50.0, epsilon = 0.001);
        assert_abs_diff_eq!(c, 10.0, epsilon = 0.001);
        assert_abs_diff_eq!(h, 180.0, epsilon = 0.001);
    }

    #[test]
    #[ignore = "This test has output which can be checked."]
    fn test_cielch_hashmap_chromaticity() {
        let gamut = CieLChGamut::new(Observer::Cie1931, CieIlluminant::D65);
        for h in 0..72 {
            print!("[");
            for l in 1..=100 {
                if let Some(c) = gamut.max_chroma(l, h) {
                    let lab = gamut.bins_to_cielab(l, c, h);
                    let [x, y] = lab.xyz().xyz().chromaticity().to_array();
                    print!("[{x:.5}, {y:.5}],");
                }
            }
            println!("],");
        }
    }
}
