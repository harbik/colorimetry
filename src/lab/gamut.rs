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

use nalgebra::DMatrix;

use crate::illuminant::CieIlluminant;
use crate::lab::{CieLCh, CieLab};
use crate::observer::Observer;
use crate::rgb::RgbSpace;
use crate::traits::Light;
use crate::xyz::XYZ;

pub struct CieLChGamut {
    illuminant: CieIlluminant,
    observer: Observer,
    max_chromas: DMatrix<u16>,
    rgb_space: RgbSpace,
    white_point: XYZ,
}

impl CieLChGamut {
    pub const H_BINS: usize = 360;
    pub const L_BINS: usize = 500;
    pub const C_BINS: usize = 500;
    const H_SCALE: f64 = Self::H_BINS as f64 / 360.0;
    const L_SCALE: f64 = Self::L_BINS as f64 / 100.0;
    const C_SCALE: f64 = Self::C_BINS as f64 / 100.0;

    pub fn new(observer: Observer, rgb_space: RgbSpace) -> Self {
        let illuminant = rgb_space.white();
        let white_point = illuminant.xyzn100(observer);
        let mut max_chromas: DMatrix<u16> = DMatrix::zeros(Self::L_BINS, Self::C_BINS);
        for hbin in 0..Self::H_BINS {
            for lbin in 0..Self::L_BINS {
                for cbin in 0..Self::C_BINS {
                    let l = Self::bin_to_l(lbin);
                    let h = Self::bin_to_h(hbin);
                    let c = Self::bin_to_c(cbin);
                    let lch = CieLCh::new([l, c, h], white_point);
                    let xyz = lch.xyz();
                    let rgb = xyz.xyz().rgb(rgb_space);
                    if rgb.is_in_gamut() {
                        max_chromas[(hbin, lbin)] = cbin as u16;
                    } else {
                        break;
                    }
                }
            }
        }
        CieLChGamut {
            illuminant,
            observer,
            max_chromas,
            white_point,
            rgb_space,
        }
    }

    pub fn rgb_space(&self) -> RgbSpace {
        self.rgb_space
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

    pub fn bins_to_cielab_static(
        whitepoint: XYZ,
        l_bin: usize,
        c_bin: usize,
        h_bin: usize,
    ) -> CieLab {
        let h = Self::bin_to_h(h_bin);
        let l = Self::bin_to_l(l_bin);
        let c = Self::bin_to_c(c_bin);
        let cielch = CieLCh::new([l, c, h], whitepoint);
        cielch.lab()
    }

    pub fn h_to_bin(value: f64) -> usize {
        (value * Self::H_SCALE).floor() as usize
    }

    pub fn bin_to_h(bin: usize) -> f64 {
        bin as f64 / Self::H_BINS as f64 * 360.0
    }

    pub fn c_to_bin(value: f64) -> usize {
        (value.clamp(0.0, 100.0) * Self::C_SCALE).floor() as usize
    }

    pub fn bin_to_c(bin: usize) -> f64 {
        bin as f64 / Self::C_BINS as f64 * 100.0
    }

    pub fn l_to_bin(value: f64) -> usize {
        (value * Self::L_SCALE).floor() as usize
    }

    pub fn bin_to_l(bin: usize) -> f64 {
        bin as f64 / Self::L_BINS as f64 * 100.0
    }

    pub fn bins_to_cielab(&self, l_bin: usize, c_bin: usize, h_bin: usize) -> CieLab {
        Self::bins_to_cielab_static(self.white_point, l_bin, c_bin, h_bin)
    }

    pub fn gamut_lch_bin(&self, lbin: usize, hbin: usize) -> CieLCh {
        let l = Self::bin_to_l(lbin);
        let h = Self::bin_to_h(hbin);
        let cbin = self.max_chromas[(hbin, lbin)] as usize;
        let c = Self::bin_to_c(cbin);
        CieLCh::new([l, c, h], self.white_point)
    }

    #[allow(unused_variables)]
    pub fn gamut_lch(&self, l: f64, h: f64) -> Result<CieLCh, crate::Error> {
        if !(0.0..100.0).contains(&l) {
            return Err(crate::Error::InvalidLightness(l));
        }
        if !(0.0..360.0).contains(&h) {
            return Err(crate::Error::InvalidHue(h));
        }

        let lbin = Self::l_to_bin(l);
        let llow = Self::bin_to_l(lbin);
        let lhigh = Self::bin_to_l(lbin + 1);
        let lfrac = (l - llow) / (lhigh - llow);

        let hbin = Self::h_to_bin(h);
        let hlow = Self::bin_to_h(hbin);
        let hhigh = Self::bin_to_h(hbin + 1);
        let hfrac = (h - hlow) / (hhigh - hlow);

        let c00 = Self::bin_to_c(self.max_chromas[(hbin, lbin)] as usize);
        let c10 = Self::bin_to_c(self.max_chromas[(hbin + 1, lbin)] as usize);
        let c01 = Self::bin_to_c(self.max_chromas[(hbin, lbin + 1)] as usize);
        let c11 = Self::bin_to_c(self.max_chromas[(hbin + 1, lbin + 1)] as usize);

        let c = (1.0 - lfrac) * ((1.0 - hfrac) * c00 + hfrac * c10)
            + lfrac * ((1.0 - hfrac) * c01 + hfrac * c11);

        Ok(CieLCh::new([l, c, h], self.white_point))
    }

    pub fn colors(&self) -> &DMatrix<u16> {
        &self.max_chromas
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_bin_conversions() {
        // Test L* conversions
        assert_eq!(
            CieLChGamut::l_to_bin(50.0),
            50 * CieLChGamut::L_SCALE as usize
        );
        assert_abs_diff_eq!(
            CieLChGamut::bin_to_l(50),
            50.0 / CieLChGamut::L_SCALE,
            epsilon = 0.001
        );
        assert_eq!(
            CieLChGamut::l_to_bin(100.0),
            100 * CieLChGamut::L_SCALE as usize
        );
        assert_abs_diff_eq!(
            CieLChGamut::bin_to_l(100),
            100.0 / CieLChGamut::L_SCALE,
            epsilon = 0.001
        );

        // Test h conversions
        assert_eq!(
            CieLChGamut::h_to_bin(180.0),
            180 * CieLChGamut::H_SCALE as usize
        );
        assert_abs_diff_eq!(
            CieLChGamut::bin_to_h(180),
            180.0 / CieLChGamut::H_SCALE,
            epsilon = 0.001
        );

        // Test C* conversions
        assert_eq!(
            CieLChGamut::c_to_bin(50.0),
            50 * CieLChGamut::C_SCALE as usize
        );
        assert_abs_diff_eq!(
            CieLChGamut::bin_to_c(50 * CieLChGamut::C_SCALE as usize),
            50.0,
            epsilon = 0.001
        );
    }

    #[test]
    fn test_white_point() {
        let gamut = CieLChGamut::new(Observer::Cie1931, RgbSpace::SRGB);
        let wp = gamut.white_point();
        assert_abs_diff_eq!(wp.x(), 95.0422, epsilon = 0.001);
        assert_abs_diff_eq!(wp.y(), 100.0, epsilon = 0.001);
        assert_abs_diff_eq!(wp.z(), 108.8610, epsilon = 0.001);
    }
}
