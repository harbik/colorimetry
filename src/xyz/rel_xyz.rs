use approx::AbsDiffEq;
use nalgebra::Vector3;

use super::XYZ;

/// # Related Tristimulus Values
/// 
/// Tristimulus Values for a given sample and reference white,
/// used to represent related colors as used in various color
/// models. Typically the reference white is normalized to have
/// an Y-value of 100 
#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct RelXYZ {
    xyz: Vector3<f64>,
    white_point: XYZ,
}

impl RelXYZ {
    pub fn new(xyz: [f64;3], white_point: XYZ) -> Self {
        RelXYZ { xyz: xyz.into(), white_point }
    }

    pub fn from_xyz(sample: XYZ, refwhite: XYZ) -> Result<Self, crate::Error> {
        if sample.observer != refwhite.observer {
            Err(crate::Error::RequireSameObserver)
        } else {
            Ok(RelXYZ { xyz: sample.xyz, white_point: refwhite })
        }
    }

    pub fn with_d65(xyz: XYZ) -> Self {
        let white_point = xyz.observer.xyz_d65();
        RelXYZ {
            xyz: xyz.xyz,
            white_point,    
        }
    }

    pub fn with_d50(xyz: XYZ) -> Self {
        let white_point = xyz.observer.xyz_d50();
        RelXYZ {
            xyz: xyz.xyz,
            white_point,    
        }
    }

    pub fn xyz(&self) -> XYZ {
        XYZ::from_vecs(self.xyz, self.white_point.observer)
    }

    pub fn white_point(&self) -> XYZ {
        self.white_point
    }

    pub fn values(&self) -> [[f64; 3]; 2] {
        [
            self.xyz.into(),
            self.white_point.xyz.into(),
        ]
    }
}

impl AbsDiffEq for RelXYZ {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let xyz_eq = self.xyz.abs_diff_eq(&other.xyz, epsilon);
        let xyzn_eq = self.white_point.abs_diff_eq(&other.white_point, epsilon);
        xyz_eq && xyzn_eq
    }
}