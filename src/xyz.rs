use nalgebra::RowVector3;
use wasm_bindgen::{prelude::wasm_bindgen, JsValue};
use crate::{obs::ObsId, Spectrum};

#[wasm_bindgen]
#[derive(Clone, Copy)]
/// A set of CIE XYZ Tristimulus values, associated with a Standard Observer.
/// They are calculated using the spectrum of a Stimulus, such as beam of light
/// reflected from a color patch, or emitted from a pixel of a display.
/// XYZ values are not often used directly, but form the basis for many colorimetric models,
/// such as CIELAB and CIECAM.
pub struct XYZ {
     pub(crate) obs: ObsId,
     pub(crate) data:  RowVector3<f64>
}


impl XYZ {

    /// Define a set of XYZ-values directly, using an identifier for its
    /// associated observer, such as `Cie::Std1931` or `Cie::Std2015`.
    pub fn new(obs: ObsId, x: f64, y: f64, z: f64) -> Self {
        let data = RowVector3::new(x,y,z);
        Self { obs: obs, data}
    }
    


    /// XYZ Tristimulus values in a [f64;3] form [X, Y, Z]
    /// ```
    /// use crate::cie::{Spectrum, CIE1931};
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = CIE1931.xyz(&Spectrum::d65()).set_illuminance(100.0);
    /// let [x, y, z] = d65_xyz.xyz();
    /// // Calculated Spreadsheet Values from CIE Datasets, over a range from 380 to 780nm
    /// assert_ulps_eq!(x, 95.042_267, epsilon = 1E-6);
    /// assert_ulps_eq!(y, 100.0);
    /// assert_ulps_eq!(z, 108.861_036, epsilon = 1E-6);
    /// ```
    pub fn xyz(&self) -> [f64; 3] {
        self.data.as_ref().clone()
    }

    pub fn set_illuminance(mut self, illuminance: f64) -> Self {
        let [x, y, z] = self.xyz();
        let s = illuminance/y;
        self.data.iter_mut().for_each(|v|*v = *v * s);
        self
    }
    
    /// Calulate Luminous value, and two dimensional (x,y) chromaticity
    /// coordinates as an array [L,x,y].
    /// ```
    /// use crate::cie::{Spectrum, CIE1931};
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = CIE1931.xyz(&Spectrum::d65());
    /// let [l,x,y] = d65_xyz.lxy();
    /// assert_ulps_eq!(x, 0.312_738, epsilon = 1E-6);
    /// assert_ulps_eq!(y, 0.329_052, epsilon = 1E-6);
    /// ```
    pub fn lxy(&self) -> [f64; 3] {
        let [x,y,z] = self.xyz();
        let s = x + y + z;
        [y, x/s, y/s]
    }

    /// CIE 1960 UCS Color Space uv coordinates *Deprecated* by the CIE, but
    /// still used for CCT calculation
    pub fn uv60(&self) -> [f64; 2] {
        let &[x, y, z] = self.data.as_ref();
        let den = x + 15.0 * y + 3.0 * z;
        [4.0 * x / den, 6.0 * y / den]
    }

    /// The CIE 1964 (U*, V*, W*) color space, also known as CIEUVW, based on
    /// the CIE 1960 UCS.
    /// Still used in Color Rendering Index Calculation
    pub fn uvw64(&self, xyz_ref: XYZ) -> [f64; 3] {
        let yy = self.data.y;
        let [ur, vr] = xyz_ref.uv60();
        let [u, v] = self.uv60();
        let ww = 25.0 * yy.powf(1.0 / 3.0) - 17.0;
        let uu = 13.0 * ww * (u - ur);
        let vv = 13.0 * ww * (v - vr);
        [uu, vv, ww]
    }

    // CIE 1976 CIELUV space, with (u',v') coordinates
    pub fn uvp(&self) -> [f64;2] {
        let &[x, y, z] = self.data.as_ref();
        let den = x + 15.0 * y + 3.0 * z;
        [4.0 * x / den, 9.0 * y / den]
    }

    pub fn cct(&self) -> [f64; 2] {
        todo!()

    }


}

impl Default for XYZ {
    fn default() -> Self {
        Self { obs: Default::default(), data: Default::default() }
    }
}

#[wasm_bindgen]
impl XYZ {

    /// XYZ Tristimuls Values JavaScript Constructor
    /// 
    /// Accepts as arguments 
    /// - x and y chromaticity coordinates only , using the "Cie::Std1931" observer as default
    /// - x and y chromaticity coordinates, and standard observer ID as 3rd argument
    /// - X, Y, and Z tristimulus values, using the "Cie::Std1931" observer as default
    /// - X, Y, and Z tristimulus values, and standard Observer ID as 4th argument

    #[wasm_bindgen(constructor)]
    pub fn new_js(x: f64, y:f64, z:JsValue, obs: JsValue) -> Self {
        todo!()
    }
}