use approx::ulps_eq;
use nalgebra::{RowVector3, Vector3};

use crate::{error::CmtError, prelude::Observer, xyz::XYZ};
use strum_macros::Display;
use wasm_bindgen::prelude::wasm_bindgen;

#[wasm_bindgen]
#[derive(Debug, Clone, Copy)]
pub struct CieLab {
    pub(crate) observer: Observer,
    pub(crate) lab: Vector3<f64>,
    pub(crate) xyzn: Vector3<f64>, // Reference white tristimulus value
}

impl CieLab {
    pub fn new(xyz: XYZ, xyzn: XYZ) -> Result<CieLab, CmtError> {
        if xyz.observer != xyzn.observer {
            Err(CmtError::RequireSameObserver)
        } else {
            Ok(CieLab {
                observer: xyz.observer,
                lab: lab(xyz.xyz, xyzn.xyz),
                xyzn: xyzn.xyz,
            })
        }
    }

    pub fn delta_e(&self, other: &Self) -> Result<f64, CmtError> {
        if ulps_eq!(self.xyzn, other.xyzn) {
            let &[l1, a1, b1] = self.lab.as_ref();
            let &[l2, a2, b2] = other.lab.as_ref();
            Ok(((l2 - l1).powi(2) + (a2 - a1).powi(2) + (b2 - b1).powi(2)).sqrt())
        } else {
            Err(CmtError::RequiresSameIlluminant)
        }
    }

    // Returns the CIE L*a*b* values as an array.
    pub fn values(&self) -> [f64; 3] {
        *self.lab.as_ref()
    }
}

impl AsRef<[f64; 3]> for CieLab {
    fn as_ref(&self) -> &[f64; 3] {
        self.lab.as_ref()
    }
}

const DELTA: f64 = 24f64 / 116f64;
const DELTA_POW2: f64 = DELTA * DELTA;
const DELTA_POW3: f64 = DELTA_POW2 * DELTA;
const LABPOW: f64 = 1f64 / 3f64;
const LABC1: f64 = 1f64 / (3f64 * DELTA_POW2);
const LABC2: f64 = 4f64 / 29f64;

fn lab_f(t: f64) -> f64 {
    if t > DELTA_POW3 {
        t.powf(LABPOW)
    } else {
        LABC1 * t + LABC2
    }
}

fn lab(xyz: Vector3<f64>, xyzn: Vector3<f64>) -> Vector3<f64> {
    let &[x, y, z] = xyz.as_ref();
    let &[xn, yn, zn] = xyzn.as_ref();
    Vector3::new(
        116f64 * lab_f(y / yn) - 16f64,
        500f64 * (lab_f(x / xn) - lab_f(y / yn)),
        200f64 * (lab_f(y / yn) - lab_f(z / zn)),
    )
}
