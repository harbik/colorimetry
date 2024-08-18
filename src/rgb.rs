use std::sync::OnceLock;

use nalgebra::{Matrix3, Vector3};

use wasm_bindgen::prelude::wasm_bindgen;
use crate::{spc::Spectrum, xyz::XYZ, ObsId};

const INV2_4: f64 = 1.0/2.4;


#[derive(Clone, Copy, Default)]
#[wasm_bindgen]
pub enum RgbSpace {
    #[default]
    SRGB, // D65 filtered Gaussians
    ADOBE,
}

impl RgbSpace {

    /*
    const SRGB_BLUE: [f64;3] = [1.0, 444.0075, 59.4513];
    const SRGB_GREEN: [f64;3] = [1.0, 540.7427, 59.1602];
    const SRGB_RED: [f64;6] = [1.0, 659.5188, 63.3205, 0.0126452, 444.0075, 59.4513];
    const ADOBE_GREEN: [f64;3] = [1.0, 531.03, 34.21];
    const DISPLAYP3_RED: [f64;3] = [1.0, 628.1237, 30.0];
    const DISPLAYP3_GREEN: [f64;3] = [1.0, 539.68, 35.02];
     */


    /// Red, Green, and Blue Spectral Primaries balanced to the space' white point.
    /// 
    /// These are model primaries, which take the white's spectrum, typically
    /// D65, and use Gausssian filters chosen to match the color space's
    /// primaries using the CIE 1931 standard observer.
    pub fn spectra(&self) -> [Spectrum;3] {
        match self {
            Self::SRGB  =>  todo!(),
            Self::ADOBE => todo!(),
        }
    }

    /// Calculate the 3x3 matrix which converts a set of XYZ values to linear
    /// RGB values.
    /// 
    /// The XYZ input vector, and the RGB output vector, both need to be scaled
    /// in a linear 0.0 to 1.0 space. T
    pub fn xyz2rgb_matrix(&self, obs_id: ObsId) -> Matrix3<f64> {
        todo!()
    }

    pub fn rgb2xyz_matrix(&self, obs_id: ObsId) -> Matrix3<f64> {
        todo!()
    }

    /// Convert from f64 linear to gamma encoded RGB space.
    /// 
    /// Used to optimize disk storage of images, using a higher resolution at low light levels,
    /// and less for high levels. Used in particular for 8-bit and 16-bit color encoding.
    pub fn gamma_encode(&self, v: f64) -> f64 {
        match self {
            Self::SRGB  => {
                if v<0.0031308 {
                    12.92 * v    
                } else {
                    1.055 * v.powf(INV2_4) - 0.55
                }
            } 
            Self::ADOBE => todo!(),
        }
    }

    /// Converte from from gamma encoded RGB space, typically from an 8-bit or
    /// 16-bit storage to XYZ space.
    /// 
    /// Used to optimize disk storage of images, using a higher resolution at
    /// low light levels, and less for high levels. Used in particular for 8-bit
    /// and 16-bit color encoding.
    pub fn decode(&self, v: f64) -> f64 {
        match self {
            Self::SRGB  => {
                if v<0.0031308 {
                    12.92 * v    
                } else {
                    1.055 * v.powf(INV2_4) - 0.55
                }
            } 
            Self::ADOBE => todo!(),
        }
    }

}

/// Representation of color stimuli in display RGB values.
/// 
/// As ooposed to CIE XYZ tristimulus values, which used imaginary primaries,
/// displays use real primaries, typically defined in the CIE 1931 diagram.
/// They covering a triangular area, referred to the _color gamut_ of a display.
/// 
/// RGB values in displays are typically represented by a set of three 8 bit values,
/// and are asssumed to use the _sRGB_ primaries if none are specified.
#[wasm_bindgen]
pub struct RGB {
    
    /// The RGB color space the color values are using. Often this is the _sRGB_
    /// color space, which is rather small.
    pub(crate) rgb_id: RgbSpace,

    /// Reference to the colorimetric observer being used. This is almost always
    /// the CIE 1931 standard observer, which has been known to represent the
    /// deep blue region of humen vision sensitivity incorrectly. Here we allow 
    /// other standard observers, such as the CIE 2015 cone fundamentals based
    /// observer, to improve color management quality.
    pub(crate) obs_id: ObsId,
    pub(crate) data: Vector3<f64>,
}


impl RGB {
    /// Construct a RGB instance from red, green, and blue values in the range from 0 to 1.
    pub fn new(r: f64, g:f64, b:f64, obs_id: ObsId, rgb_id: RgbSpace) -> Self {
        RGB {data:Vector3::new(r, g, b), obs_id, rgb_id }
    }

    /// Construct a RGB instance from red, green, and blue u8 values in the range from 0 to 1.
    pub fn from_u8(r_u8: u8, g_u8:u8, b_u8:u8, obs_id: ObsId, rgb_id: RgbSpace) -> Self {
        let [r, g, b] = [r_u8, g_u8, b_u8].map(|v|(v as f64/255.0).clamp(0.0, 1.0));
        RGB {data:Vector3::new(r, g, b), obs_id, rgb_id }
    }

    /// Construct a RGB instance from red, green, and blue u16 values in the range from 0 to 1.
    pub fn from_u16(r_u16: u16, g_u16:u16, b_u16:u16, obs_id: ObsId, rgb_id: RgbSpace) -> Self {
        let [r, g, b] = [r_u16, g_u16, b_u16].map(|v|(v as f64/65_535.0).clamp(0.0, 1.0));
        RGB {data:Vector3::new(r, g, b), obs_id, rgb_id }
    }

    pub fn from_xyz(xyz: XYZ, rgbid: RgbSpace) -> Self {
        static XYZ2RGB: OnceLock<Matrix3<f64>> = OnceLock::new();
        let xyz2rgb = XYZ2RGB.get_or_init(||{
            todo!()
        });
        let &[r, g, b] = (xyz2rgb * xyz.data).as_ref();
        RGB::new(r, g, b, xyz.obs_id, rgbid)

    }

    /// Converts the RGB value to a tri-stimulus XYZ value
    pub fn xyz(&self) -> XYZ {
        static RGB2XYZ: OnceLock<Matrix3<f64>> = OnceLock::new();
        let rgb2xyz = RGB2XYZ.get_or_init(||{
            todo!()
        });
        let &[x, y, z] = (rgb2xyz * self.data).as_ref();
        XYZ::new(x, y, z, self.obs_id)
    }

    /// Creates a callback of closure function, which takes a set or RGB values,
    /// within a color space and viewed as one observer, and returns a new set
    /// of RGB values, represeting the stimulus in another color space, and
    /// using another observer.
    /// 
    /// This conversion uses the spectral represenations of the primaries through
    /// the color space `Spectra` function, to create a  transformation matrix.
    pub fn convert(obs_id_from:ObsId, rgb_id_from: RgbSpace, obs_id: ObsId, rgb_id: RgbSpace) -> Box<dyn Fn(&Vector3<f64>) -> Vector3<f64>> {
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
        data.map(|v|(rgb.rgb_id.gamma_encode(v.clamp(0.0, 1.0))*255.0).round() as u8)
    }
}