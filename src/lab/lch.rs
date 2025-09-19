// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

use approx::{ulps_eq, AbsDiffEq, UlpsEq};
use nalgebra::Vector3;

use crate::{
    lab::CieLab,
    observer::Observer,
    rgb::{RgbSpace, WideRgb},
    xyz::{RelXYZ, XYZ},
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CieLCh {
    lch: Vector3<f64>,
    white_point: XYZ,
}

impl CieLCh {
    pub fn new(lch: [f64; 3], xyzn: XYZ) -> Self {
        Self {
            lch: lch.into(),
            white_point: xyzn,
        }
    }

    pub fn l(&self) -> f64 {
        self.lch[0]
    }

    pub fn c(&self) -> f64 {
        self.lch[1]
    }

    pub fn h(&self) -> f64 {
        self.lch[2]
    }

    pub fn white_point(&self) -> XYZ {
        self.white_point
    }

    pub fn observer(&self) -> Observer {
        self.white_point.observer()
    }

    pub fn to_array(&self) -> [f64; 3] {
        self.lch.into()
    }

    pub fn from_xyz(rxyz: RelXYZ) -> Self {
        let lab = CieLab::from_rxyz(rxyz);
        Self::from_lab(lab)
    }

    pub fn from_lab(lab: CieLab) -> Self {
        let &[l, a, b] = lab.as_ref();
        let c = (a * a + b * b).sqrt();
        // if a and b are both zero, atan2 will return 0.0, which is correct for hue
        let h = b.atan2(a).to_degrees().rem_euclid(360.0); // Convert to degrees and normalize
        Self {
            lch: Vector3::new(l, c, h),
            white_point: lab.white_point(),
        }
    }

    /// This method converts the CIE LCh color representation to CIE Lab.
    pub fn lab(&self) -> CieLab {
        let [l, c, h] = self.lch.into();
        let h_rad = h.to_radians();
        let a = c * h_rad.cos();
        let b = c * h_rad.sin();
        CieLab::new([l, a, b], self.white_point)
    }

    /// This method converts the CIE LCh color representation to related XYZ values.
    pub fn rxyz(&self) -> RelXYZ {
        self.lab().rxyz()
    }

    /// This method converts the CIE LCh color representation to related XYZ values.
    pub fn xyz(&self) -> XYZ {
        self.rxyz().xyz()
    }

    pub fn rgb(&self, rgb_space: RgbSpace) -> WideRgb {
        // Convert CIE LCh to CIE Lab, then to XYZ, and finally to RGB
        self.xyz().rgb(rgb_space)
    }

    /// Check the validity of the CIE LCh color representation.
    /// Returns `true` if the values are within valid ranges and the conversion back its Relative XYZ values is consistent.
    pub fn is_valid(&self) -> bool {
        // Values should be not negative, L should be in [0, 100], and C should be non-negative
        if !self.lch.iter().all(|&v| v.is_finite() && v >= 0.0)
            || self.lch[0] > 100.0
            || self.lch[2] > 360.0
        {
            return false;
        }

        // Convert to XYZ and check if it is valid
        let xyz = self.rxyz();
        if !xyz.is_valid() {
            return false;
        }

        // Convert back to CIE LCh and check if it matches the original
        let lch_back = Self::from_xyz(xyz);
        ulps_eq!(self, &lch_back)
    }

    pub fn ciede(&self, other: &Self) -> Result<f64, crate::Error> {
        // Calculate the CIEDE2000 color difference between two CIE LCh colors
        let lab1 = self.lab();
        let lab2 = other.lab();
        lab1.ciede(&lab2)
    }

    pub fn ciede2000(&self, other: &Self) -> Result<f64, crate::Error> {
        // Calculate the CIEDE2000 color difference between two CIE LCh colors
        let lab1 = self.lab();
        let lab2 = other.lab();
        lab1.ciede2000(&lab2)
    }
}

/// Implment the AbsDiffEq and UlpsEq traits for CieLCh to allow for approximate equality checks
/// using absolute difference and ULPs (Units in the Last Place) respectively.
impl AbsDiffEq for CieLCh {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-6
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.lch.abs_diff_eq(&other.lch, epsilon)
            && self.white_point.abs_diff_eq(&other.white_point, epsilon)
    }
}

impl UlpsEq for CieLCh {
    fn default_max_ulps() -> u32 {
        10
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.lch.ulps_eq(&other.lch, epsilon, max_ulps)
            && self
                .white_point
                .ulps_eq(&other.white_point, epsilon, max_ulps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::illuminant::CieIlluminant;
    use crate::observer::Observer::Cie1931;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_lch_conversion() {
        let lch = CieLCh::new([50.0, 20.0, 30.0], Cie1931.xyz_d65());
        let xyz = lch.rxyz();

        assert_abs_diff_eq!(lch.l(), 50.0, epsilon = 1e-6);
        assert_abs_diff_eq!(lch.c(), 20.0, epsilon = 1e-6);
        assert_abs_diff_eq!(lch.h(), 30.0, epsilon = 1e-6);
        assert!(lch.is_valid());
        assert!(xyz.is_valid());
    }

    #[test]
    fn test_roundtrip_monochromes() {
        let monos = Cie1931.monochromes(CieIlluminant::D65);
        monos.into_iter().for_each(|(_l, xyz)| {
            let lch = CieLCh::from_xyz(xyz);
            let xyz = lch.rxyz();
            let lch_back = CieLCh::from_xyz(xyz);
            assert_abs_diff_eq!(lch, lch_back, epsilon = 1e-10);
        });
    }
}
