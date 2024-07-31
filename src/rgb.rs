use nalgebra::{Matrix3, RowVector3};

use wasm_bindgen::prelude::wasm_bindgen;
use crate::spc::Spectrum;


#[wasm_bindgen]
pub enum RGBSpace {
    SRGB, //  0..255
    SRGBFloat, // 0.0..1.0
    SRGB16, // 0..64K
//    AdobeRGB,
}

impl RGBSpace {

    pub fn scale(&self, v: f64) -> f64 {
        match self {
            RGBSpace::SRGB => v.clamp(0.0, 255.0)/255.0,
            RGBSpace::SRGBFloat => v,
            RGBSpace::SRGB16 => v.clamp(0.0, 65_535.0)/65_535.0,
        }
    }

    pub fn xyz2rgb(&self) -> Matrix3<f64> {
        todo!()
    }

    pub fn rgb2xyz(&self) -> Matrix3<f64> {
        todo!()
    }

    /*
    const SRGB_BLUE: [f64;3] = [1.0, 444.0075, 59.4513];
    const SRGB_GREEN: [f64;3] = [1.0, 540.7427, 59.1602];
    const SRGB_RED: [f64;6] = [1.0, 659.5188, 63.3205, 0.0126452, 444.0075, 59.4513];
    const ADOBE_GREEN: [f64;3] = [1.0, 531.03, 34.21];
    const DISPLAYP3_RED: [f64;3] = [1.0, 628.1237, 30.0];
    const DISPLAYP3_GREEN: [f64;3] = [1.0, 539.68, 35.02];
     */


    pub fn spectra(&self) -> [Spectrum;3] {
        /*
        match self {
            self::SRGB | self::SRGBFloat | self::SRGB16 => [
                Spectrum::led(center, width),
                Spectrum::led(center, width),
                Spectrum::led(center, width) + Spectrum::led(center, width),
                ],

        }
        */
        todo!()
    }
}

#[wasm_bindgen]
pub struct RGB {
    pub(crate) id: RGBSpace,
    pub(crate) data: RowVector3<f64>,
}


impl RGB {
    pub fn new(r: f64, g:f64, b:f64, id: RGBSpace) -> Self {
        RGB {data:RowVector3::new(id.scale(r), id.scale(g), id.scale(b)), id }
    }

    pub fn xyz(&self) -> Self {
        todo!()
    }

}
