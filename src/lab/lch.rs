use approx::{abs_diff_eq, AbsDiffEq};
use nalgebra::Vector3;

use crate::{
    lab::CieLab,
    xyz::{RelXYZ, XYZ},
};

#[derive(Debug, Clone, Copy)]
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

    pub fn lch(&self) -> [f64; 3] {
        self.lch.into()
    }

    pub fn l(&self) -> f64 {
        self.lch()[0]
    }

    pub fn c(&self) -> f64 {
        self.lch()[1]
    }

    pub fn h(&self) -> f64 {
        self.lch()[2]
    }

    pub fn white_point(&self) -> XYZ {
        self.white_point
    }

    pub fn values(&self) -> [f64; 3] {
        self.lch.into()
    }

    pub fn from_xyz(rxyz: RelXYZ) -> Self {
        let lab = CieLab::from_xyz(rxyz);
        Self::from_lab(lab)
    }

    pub fn from_lab(lab: CieLab) -> Self {
        let &[l, a, b] = lab.as_ref();
        let c = (a * a + b * b).sqrt();
        let h = if c < 0.001 {
            0.0 // set to 0 if chroma is negligible
        } else {
            b.atan2(a).to_degrees().rem_euclid(360.0) // Convert to degrees and normalize
        };
        Self {
            lch: Vector3::new(l, c, h),
            white_point: lab.xyzn(),
        }
    }

    pub fn lab(&self) -> CieLab {
        let [l, c, h] = self.lch();
        let h_rad = h.to_radians();
        let a = c * h_rad.cos();
        let b = c * h_rad.sin();
        CieLab::new([l, a, b], self.white_point)
    }

    pub fn xyz(&self) -> RelXYZ {
        self.lab().xyz()
    }

    pub fn is_valid(&self) -> bool {
        let xyz = self.xyz();
        if !xyz.is_valid() {
            return false;
        }
        let lch_back = Self::from_xyz(xyz);
        abs_diff_eq!(*self, lch_back, epsilon = 1e-6)
    }
}

impl PartialEq for CieLCh {
    fn eq(&self, other: &Self) -> bool {
        self.lch == other.lch && self.white_point == other.white_point
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::illuminant::CieIlluminant;
    use crate::observer::Observer::Cie1931;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_lch_conversion() {
        let lch = CieLCh::new([50.0, 20.0, 30.0], Cie1931.xyz_d65());
        let xyz = lch.xyz();

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
            let xyz = lch.xyz();
            let lch_back = CieLCh::from_xyz(xyz);
            assert_abs_diff_eq!(lch, lch_back, epsilon = 1e-10);
        });
    }
}
