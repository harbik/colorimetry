use nalgebra::RowVector3;

use strum_macros::Display;
use wasm_bindgen::prelude::wasm_bindgen;
use crate::xyz::XYZ;

#[wasm_bindgen]
#[derive(Debug, Clone, Copy)]
pub struct CieLab {
    pub(crate) data: RowVector3<f64>,
    pub(crate) xyz_white: XYZ, // contains obs
}

impl CieLab {
    pub fn new(x: f64, y:f64, z:f64, xyz_white: XYZ) -> CieLab {
        CieLab {data: lab(x, y, z, xyz_white), xyz_white }
    }
}

impl AsRef<[f64;3]> for CieLab {
    fn as_ref(&self) -> &[f64;3] {
        self.data.as_ref()
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

fn lab(x: f64, y: f64, z: f64, xyz_r: XYZ) -> RowVector3<f64> {
    let &[xr, yr, zr] = xyz_r.xyz.as_ref();
    RowVector3::new(
        116f64 * lab_f(y/yr) - 16f64, 
        500f64 * (lab_f(x/xr) - lab_f(y/yr)), 
        200f64 * (lab_f(y/yr) - lab_f(z/zr)) 
    )

}