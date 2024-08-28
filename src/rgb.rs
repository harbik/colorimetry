use std::{collections::HashMap, primitive, sync::{LazyLock, OnceLock}};

use nalgebra::{Matrix3, Vector3};

use strum_macros::EnumIter;
use wasm_bindgen::prelude::wasm_bindgen;
use crate::{gamma::GammaCurve, spc::Spectrum, xyz::XYZ, ObsId, StdIlluminant, CIE1931, D65};

const INV2_4: f64 = 1.0/2.4;

// The display P3 red coordinate is outside the CIE 1931 gamut using the CIE
// 1931 1 nanometer dataset as provided by the CIE.  To still match it, I need
// to desature it by mixing in a bit of white. It mxes in a bit of blue (1.4E-6)
// to get to the deaturated target.
const D:f64 = 0.0005620; // Desaturation ratio
const D65X:f64 = 0.312_738;
const D65Y:f64 = 0.329_052;

/// These are the Chromaticities of the various RGB colorspaces, using the CIE 1931 standaard observer.
/// This uses `LazyLock` to create a global static varibale, which gets intialized when first referenced.
pub static XY_PRIMARIES: LazyLock<HashMap<&str, ([[f64;2];3], StdIlluminant)>> = LazyLock::new(|| {
        HashMap::from([
            ("sRGB", ([[0.64, 0.33], [0.3, 0.6], [0.15, 0.06]], StdIlluminant::D65)),
            ("Adobe RGB", ([[0.64, 0.33], [0.21, 0.71], [0.15, 0.06]], StdIlluminant::D65)),
            ("Display P3", ([[(1.0-D)*0.68+D*D65X, (1.0-D)*0.32+D*D65Y], [0.265, 0.69], [0.15, 0.06]], StdIlluminant::D65)),
        ])
    }
);


#[derive(Debug, Clone, Copy, Default, EnumIter)]
#[wasm_bindgen]
/**
A Light Weight index tag, to represent an RGB space.
Used for example in the RGB value set, to identify the color space being used.  
 */
pub enum RgbSpaceId {
    #[default]
    SRGB, // D65 filtered Gaussians
    ADOBE,
    DisplayP3,
}

impl RgbSpaceId {

    /**
    Obtain reference to the RgbSpace data, and a color space name string.

    An `RgbSpace` contains the primary and white spectra, and the gamma
    curve function used for encoding and decoding RGB values into a data
    stream or image pixel.
    */
    pub fn rgb_space(&self) -> (&RgbSpace, &str) {
        match self {
            Self::SRGB  =>  (RgbSpace::srgb(),"sRGB"),
            Self::ADOBE => (RgbSpace::adobe_rgb(), "Adobe RGB"),
            Self::DisplayP3 => (RgbSpace::display_p3(), "Display P3"),
        }

    }
}


/**
Spectrally based color space, using spectral representations of the primaries and the
reference white.

Using the CIE 1931 standard observer, using a wavelength domain from 380 top 780
nanometer with 1 nanometer steps, these result in their usual chromaticity
values.  The most common _sRGB_ color space is obtained using the
`RgbSpace::srgb()` constructor. For this instance, the blue and green primaries
are direct Gaussian-filtered D65 spectra. A mixture of the blue primary and a
G1aussian-filtered red component is used for the red primary. Similar
constructors are provided for the `Adobe` and `DisplayP3` color spaces.

The benefit of spectral primaries is that color management and color profiles
can use updated Colorimetric Observers, such as the Cone-Fundamental based CIE
2015 observers, which don't have the CIE 1931 deficiencies. For example, they
can also be optimized for special observers by considering an observer's age or
health conditions.
*/
pub struct RgbSpace {
   // pub(crate) id: RgbSpaceId,
    pub(crate) primaries: [Spectrum;3],
    pub(crate) white: Spectrum,
    pub(crate) gamma: GammaCurve,
}


impl RgbSpace {
    
    /**
    Creates a new `RgbSpace` using an array of three spectra as first
    argument, representing the spectral distributions of the red, green, and
    blue primary stimuli, the white spectrum as second argument, and a gamma
    function as last arguments.  The spectral primaries or white spectrum
    can have arbitrary power, as they are normalized before using in
    calculation, in particular when using the RGB to XYZ transforms, and
    vice versa (see the `rgb2xyz(rgbid: &RgbSpaceId)` and `xyz2rgb(rgbid: &RgbSpaceId)`
    methods of `Observer`).
    */
    pub fn new(primaries: [Spectrum;3], white: Spectrum, gamma: GammaCurve) -> Self {
        Self { primaries, white, gamma }
    }


    /**
    The sRGB color space, created by HP and Microsoft in 1996.  It is the
    default color space used in an image, or a Web-page, if no color space
    is specified in a tag or in a color profile.
    In this instance, the primaries are composed of Daylight D65 filtered Gaussian filters.
    For the green and blue these are single Gaussians.
    For red primary a combination of a blue and red Gaussian is used, with
    the blue being equal to the blue primary.
    */
    pub fn srgb()-> &'static RgbSpace {
    static SRGB: OnceLock<RgbSpace> = OnceLock::new();
        /*
        Red Gaussian center wavelength and width, and blue lumen fraction
        added to match the red primary sRGB color space specifiation.
        See the `examples/primaries/main.rs` how these were calculated.
        */
        const RED: [f64;3] = [627.2101041540204, 23.38636113607498,0.006839349789397155];

        // Blue and Green Gaussian filters.
        const GREEN: [f64;2] = [541.2355906535001, 33.66683554608116];
        const BLUE: [f64;2] = [398.0273721579992, 52.55338039394701];

        SRGB.get_or_init(||{
            let primaries = gaussian_filtered_primaries(&D65, RED, GREEN, BLUE);
            let white = Spectrum::d65_illuminant();
            let gamma = GammaCurve::new(vec![2.4, 1.0/1.055, 0.055/1.055, 1.0/12.92, 0.04045]);
            Self { primaries, white, gamma}
        })
    }
    /**
    The Adobe RGB color space.

    The primaries are composed of Daylight D65 filtered Gaussian filters.
    For the green and blue these are single Gaussians.
    For red primary a combination of a blue and red Gaussian is used, with
    the blue being equal to the blue primary.
    */
    pub fn adobe_rgb()-> &'static RgbSpace {
    static ADOBE_RGB: OnceLock<RgbSpace> = OnceLock::new();
        // Red Gaussian center wavelength and width, and blue lumen fraction
        // added to match the red primary sRGB color space specifiation.
        // See the `examples/primaries/main.rs` how these were calculated.
        const RED: [f64;3] = [627.2101041540204, 23.38636113607498, 0.006839349789397155];

        // Blue and Green Gaussian filters.
        const GREEN: [f64;2] = [531.1505632157933, 20.61013919689458];
        const BLUE: [f64;2] = [398.0273721579992, 52.55338039394701];

        ADOBE_RGB.get_or_init(||{
            let primaries = gaussian_filtered_primaries(&D65, RED, GREEN, BLUE);
            let white = Spectrum::d65_illuminant();
            let gamma = GammaCurve::new(vec![563.0/256.0]);
                // See https://en.wikipedia.org/wiki/Adobe_RGB_color_space#ICC_PCS_color_image_encoding
            Self { primaries, white, gamma}
        })
    }

    /**
    The Display P3 color space.
    The primaries are composed of Daylight D65 filtered Gaussian filters.
    For the green and blue these are single Gaussians.
    For red primary a combination of a blue and red Gaussian is used, with
    the blue being equal to the blue primary.
    */
    pub fn display_p3()-> &'static RgbSpace {
    static DISPLAY_P3: OnceLock<RgbSpace> = OnceLock::new();
        // Red Gaussian center wavelength and width, and blue lumen fraction
        // added to match the red primary sRGB color space specifiation.
        // See the `examples/primaries/main.rs` how these were calculated.
        const RED: [f64;3] = [637.9658554073235, 23.274151215539906, 1.4261447449730065e-6];

        // Blue and Green Gaussian filters.
        const GREEN: [f64;2] = [539.8416064376327, 21.411199475289777];
        const BLUE: [f64;2] = [398.0273721579992, 52.55338039394701];

        DISPLAY_P3.get_or_init(||{
            let primaries = gaussian_filtered_primaries(&D65, RED, GREEN, BLUE);
            let white = Spectrum::d65_illuminant();
            let gamma = GammaCurve::new(vec![2.4, 1.0/1.055, 0.055/1.055, 1.0/12.92, 0.04045]);
            Self { primaries, white, gamma}
        })
    }
}

#[cfg(test)]
mod rgbspace_tests {
    use crate::{rgb::RgbSpace, RgbSpaceId, CIE1931, XY_PRIMARIES, Spectrum, D65};
    use approx::assert_ulps_eq;
    use strum::IntoEnumIterator;

    #[test]
    /// Check color points of the primaries, as calculated from the space's
    /// spectra, to the targets in `XY_PRIMARIES`. 
    fn srgb_test(){
        for space in RgbSpaceId::iter() {
            let (rgbspace, rgbstr) = space.rgb_space();
            for i in 0..3 {
                let xy = CIE1931.xyz(&rgbspace.primaries[i]).chromaticity();
                let xywant = XY_PRIMARIES[rgbstr].0[i];
                assert_ulps_eq!(xy.as_ref(), xywant.as_ref(), epsilon = 1E-5);

            }
        }
    }

}


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
    pub(crate) rgb_id: RgbSpaceId,

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
    pub fn new(r: f64, g:f64, b:f64, obs_id: ObsId, rgb_id: RgbSpaceId) -> Self {
        RGB {data:Vector3::new(r, g, b), obs_id, rgb_id }
    }

    /// Construct a RGB instance from red, green, and blue u8 values in the range from 0 to 1.
    pub fn from_u8(r_u8: u8, g_u8:u8, b_u8:u8, obs_id: ObsId, rgb_id: RgbSpaceId) -> Self {
        let [r, g, b] = [r_u8, g_u8, b_u8].map(|v|(v as f64/255.0).clamp(0.0, 1.0));
        RGB {data:Vector3::new(r, g, b), obs_id, rgb_id }
    }

    /// Construct a RGB instance from red, green, and blue u16 values in the range from 0 to 1.
    pub fn from_u16(r_u16: u16, g_u16:u16, b_u16:u16, obs_id: ObsId, rgb_id: RgbSpaceId) -> Self {
        let [r, g, b] = [r_u16, g_u16, b_u16].map(|v|(v as f64/65_535.0).clamp(0.0, 1.0));
        RGB {data:Vector3::new(r, g, b), obs_id, rgb_id }
    }

    pub fn from_xyz(xyz: XYZ, rgbid: RgbSpaceId) -> Self {
        static XYZ2RGB: OnceLock<Matrix3<f64>> = OnceLock::new();
        let xyz2rgb = XYZ2RGB.get_or_init(||{
            todo!()
        });
        let &[r, g, b] = (xyz2rgb * xyz.data).as_ref();
        RGB::new(r, g, b, xyz.obs_id, rgbid)

    }

    /// Converts the RGB value to a tri-stimulus XYZ value
    pub fn xyz(&self) -> XYZ {
        const YW: f64 = 100.0;
        let data = self.obs_id.observer().rgb2xyz(&self.rgb_id) * self.data;
        let s = YW/data.y;
        XYZ {
            obs_id: self.obs_id,
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
    pub fn convert(obs_id_from:ObsId, rgb_id_from: RgbSpaceId, obs_id: ObsId, rgb_id: RgbSpaceId) -> Box<dyn Fn(&Vector3<f64>) -> Vector3<f64>> {
        todo!()

    }

    /// Transform a set of RGB values, defining a stimulus for one standard
    /// observer, into a set of RGB values representing the same stimulus for
    /// different standard observer or special observer.  On initial use this
    /// function calculates a transformation matrix based on the colorimetric
    /// tristimulus values of the respective primaries.
    pub fn transform(&self, obs_from: &ObsId) -> Self {
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
        data.map(|v|(rgb.rgb_id.rgb_space().0.gamma.encode(v.clamp(0.0, 1.0))*255.0).round() as u8)
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
