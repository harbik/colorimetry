use std::sync::OnceLock;
use nalgebra::{Matrix3, Vector3};
use wasm_bindgen::prelude::wasm_bindgen;
use crate::{spc::Spectrum, xyz::XYZ, ObserverTag, RgbSpaceTag, CIE1931};


/// Representation of a color stimulus in a set of Red, Green, and Blue (RGB) values,
/// representing its relative composition using standard primaries.
/// 
/// RGB values are commonly used in digital images, with the relative intensity
/// of the primaries defined as three 8-bit values, with range from 0 to 255.
/// As ooposed to CIE XYZ tristimulus values, which used imaginary primaries,
/// displays use real primaries, typically defined in the CIE 1931 diagram.
/// They cover a triangular area, referred to the _color gamut_ of a display.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy)]
pub struct RGB {
    
    /// The RGB color space the color values are using. Often this is the _sRGB_
    /// color space, which is rather small.
    pub(crate) space: RgbSpaceTag,

    /// Reference to the colorimetric observer being used. This is almost always
    /// the CIE 1931 standard observer, which has been known to represent the
    /// deep blue region of humen vision sensitivity incorrectly. Here we allow 
    /// other standard observers, such as the CIE 2015 cone fundamentals based
    /// observer, to improve color management quality.
    pub(crate) obs: ObserverTag,
    pub(crate) data: Vector3<f64>,
}


impl RGB {
    /// Construct a RGB instance from red, green, and blue values in the range from 0 to 1.
    pub fn new(r: f64, g:f64, b:f64, obs: ObserverTag, space: RgbSpaceTag) -> Self {
        RGB {data:Vector3::new(r, g, b), obs, space }
    }

    /// Construct a RGB instance from red, green, and blue u8 values in the range from 0 to 1.
    pub fn from_u8(r_u8: u8, g_u8:u8, b_u8:u8, obs: ObserverTag, space: RgbSpaceTag) -> Self {
        let [r, g, b] = [r_u8, g_u8, b_u8].map(|v|(v as f64/255.0).clamp(0.0, 1.0));
        RGB {data:Vector3::new(r, g, b), obs, space }
    }

    /// Construct a RGB instance from red, green, and blue u16 values in the range from 0 to 1.
    pub fn from_u16(r_u16: u16, g_u16:u16, b_u16:u16, obs: ObserverTag, space: RgbSpaceTag) -> Self {
        let [r, g, b] = [r_u16, g_u16, b_u16].map(|v|(v as f64/65_535.0).clamp(0.0, 1.0));
        RGB {data:Vector3::new(r, g, b), obs, space }
    }

    pub fn from_xyz(xyz: XYZ, space: RgbSpaceTag) -> Self {
        static XYZ2RGB: OnceLock<Matrix3<f64>> = OnceLock::new();
        let xyz2rgb = XYZ2RGB.get_or_init(||{
            todo!()
        });
        let &[r, g, b] = (xyz2rgb * xyz.data).as_ref();
        RGB::new(r, g, b, xyz.obs_id, space)

    }

    /// Converts the RGB value to a tri-stimulus XYZ value
    pub fn xyz(&self) -> XYZ {
        const YW: f64 = 100.0;
        let data = self.obs.observer().rgb2xyz(&self.space) * self.data;
        let s = YW/data.y;
        XYZ {
            obs_id: self.obs,
            data: data.map(|v|v*s),
            yw: Some(YW)
        }
    }

    /// Creates a callback of closure function, which takes a set or RGB values,
    /// within a color space and viewed as one observer, and returns a new set
    /// of RGB values, represeting the stimulus in another color space, and
    /// using another observer.
    /// 
    /// This conversion uses the spectral represenations of the primaries through
    /// the color space `Spectra` function, to create a  transformation matrix.
    pub fn convert(obs_from:ObserverTag, space_from: RgbSpaceTag, obs: ObserverTag, space: RgbSpaceTag) -> Box<dyn Fn(&Vector3<f64>) -> Vector3<f64>> {
        todo!()

    }

    /// Transform a set of RGB values, defining a stimulus for one standard
    /// observer, into a set of RGB values representing the same stimulus for
    /// different standard observer or special observer.  On initial use this
    /// function calculates a transformation matrix based on the colorimetric
    /// tristimulus values of the respective primaries.
    pub fn transform(&self, obs_from: &ObserverTag) -> Self {
        todo!() 
    }


}

impl AsRef<Vector3<f64>> for RGB {
    fn as_ref(&self) -> &Vector3<f64> {
        &self.data
    }
}

/// Clamped RGB values as a u8 array. Uses gamma function.
impl From<RGB> for [u8;3] {
    fn from(rgb: RGB) -> Self {
        let data: &[f64;3] =rgb.data.as_ref();
        data.map(|v|(rgb.space.rgb_space().0.gamma.encode(v.clamp(0.0, 1.0))*255.0).round() as u8)
    }
}

pub fn gaussian_filtered_primaries(white: &Spectrum, red: [f64;3], green: [f64;2], blue: [f64;2]) -> [Spectrum; 3] {
    let [rc, rw, f] = red;
    let [gc, gw]=  green;
    let [bc, bw]=  blue;
    [
        Spectrum::gaussian_filter(bc, bw).mul(white).set_illuminance(&CIE1931, 100.0).mul_f64(f) +
        Spectrum::gaussian_filter(rc, rw).mul(white).set_illuminance(&CIE1931, 100.0).mul_f64(1.0-f),
        Spectrum::gaussian_filter(gc, gw).mul(white).set_illuminance(&CIE1931, 100.0),
        Spectrum::gaussian_filter(bc, bw).mul(white).set_illuminance(&CIE1931, 100.0),
    ]
}
