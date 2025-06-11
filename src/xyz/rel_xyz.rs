use approx::AbsDiffEq;

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
    related: XYZ,
    white_point: XYZ,
}

impl RelXYZ {
    pub fn new(sample: XYZ, refwhite: XYZ) -> Result<Self, crate::Error> {
        if sample.observer != refwhite.observer {
            Err(crate::Error::RequireSameObserver)
        } else {
            Ok(RelXYZ { related: sample, white_point: refwhite })
        }
    }

    pub fn with_d65(related: XYZ) -> Self {
        let white_point = related.observer.xyz_d65();
        RelXYZ {
            related,
            white_point,    
        }
    }

    pub fn with_d50(related: XYZ) -> Self {
        let white_point = related.observer.xyz_d50();
        RelXYZ {
            related,
            white_point,    
        }
    }

    pub fn xyz(&self) -> XYZ {
        self.related
    }

    pub fn white_point(&self) -> XYZ {
        self.white_point
    }

    pub fn values(&self) -> [[f64; 3];2] {
        [self.related.values(), self.white_point.values()]
    }
}

impl AbsDiffEq for RelXYZ {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let xyz_eq = self.related.abs_diff_eq(&other.related, epsilon);
        let xyzn_eq = self.white_point.abs_diff_eq(&other.white_point, epsilon);
        let obs_eq = self.related.observer == other.related.observer;
        xyz_eq && xyzn_eq && obs_eq
    }
}