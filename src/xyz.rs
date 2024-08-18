use nalgebra::Vector3;
use crate::obs::ObsId;
use wasm_bindgen::prelude::wasm_bindgen; 


#[wasm_bindgen]
#[derive(Clone, Copy)]
/// A set of CIE XYZ Tristimulus values, associated with a Standard Observer.
/// They are calculated using the spectrum of a Stimulus, such as beam of light
/// reflected from a color patch, or emitted from a pixel of a display.
/// XYZ values are not often used directly, but form the basis for many colorimetric models,
/// such as CIELAB and CIECAM.
pub struct XYZ {
     pub(crate) obs_id: ObsId,
     pub(crate) data:  Vector3<f64>
}


impl XYZ {

    /// Define a set of XYZ-values directly, using an identifier for its
    /// associated observer, such as `Cie::Std1931` or `Cie::Std2015`.
    pub fn new(x: f64, y: f64, z: f64, obs_id: ObsId) -> Self {
        let data = Vector3::new(x,y,z);
        Self { obs_id, data}
    }
    


    /// XYZ Tristimulus values in a [f64;3] form [X, Y, Z]
    /// ```
    /// use crate::colorimetry::{Spectrum, CIE1931};
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = CIE1931.xyz(&Spectrum::d65_illuminant()).set_illuminance(100.0);
    /// let [x, y, z] = d65_xyz.xyz();
    /// // Calculated Spreadsheet Values from CIE Datasets, over a range from 380 to 780nm
    /// assert_ulps_eq!(x, 95.042_267, epsilon = 1E-6);
    /// assert_ulps_eq!(y, 100.0);
    /// assert_ulps_eq!(z, 108.861_036, epsilon = 1E-6);
    /// ```
    pub fn values(&self) -> [f64; 3] {
        self.data.as_ref().clone()
    }

    pub fn set_illuminance(mut self, illuminance: f64) -> Self {
        let [_x, y, _z] = self.values();
        let s = illuminance/y;
        self.data.iter_mut().for_each(|v|*v = *v * s);
        self
    }
    
    /// The chromaticity coordinates as an array [x,y].
    /// ```
    /// use crate::colorimetry::{Spectrum, CIE1931};
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = CIE1931.xyz(&Spectrum::d65_illuminant());
    /// let [l,x,y] = d65_xyz.lxy();
    /// assert_ulps_eq!(x, 0.312_738, epsilon = 1E-6);
    /// assert_ulps_eq!(y, 0.329_052, epsilon = 1E-6);
    /// ```
    pub fn chromaticity(&self) -> [f64; 2] {
        let [x,y,z] = self.values();
        let s = x + y + z;
        [x/s, y/s]
    }

    /// Luminous value Y of the tristimulus values.
    /// 
    /// This value is associated with different types of photometric quantities -
    /// see Wikipedia's article on [Luminous Intensity](https://en.wikipedia.org/wiki/Luminous_intensity).
    /// As a stimuilus, this value has a unit of candela per square meters,
    /// as an illuminant lux, or lumen per square meter.
    pub fn luminous_value(&self) -> f64 {
        self.data.y
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
        Self { obs_id: Default::default(), data: Default::default() }
    }
}

// JS-WASM Interface code
#[cfg(target_arch="wasm32")]
#[wasm_bindgen]
impl XYZ {

    /**
    Create an XYZ Tristimuls Values object.

    Accepts as arguments 

    - x and y chromaticity coordinates only , using the "Cie::Std1931" observer as default
    - x and y chromaticity coordinates, and standard observer ID as 3rd argument
    - X, Y, and Z tristimulus values, using the "Cie::Std1931" observer as default
    - X, Y, and Z tristimulus values, and a standard Observer ID as 4th argument

    When only x and y chromaticity coordinates are specified, the luminous
    value is set to 100.0 candela per square meter.

    ```javascript, ignore
    // Create a new XYZ object using D65 CIE 1931 chromaticity coordinates
    const xyz = new cmt.XYZ(0.31272, 0.32903);
    
    // Get and check the corresponding tristimulus values, with a luminous value
    // of 100.0
    const [x, y, z] = xyz.values();
    assert.assertAlmostEquals(x, 95.047, 5E-3); // D65 wikipedia
    assert.assertAlmostEquals(y, 100.0);
    assert.assertAlmostEquals(z, 108.883, 5E-3);

    // and get back the orgiinal chromaticity coordinates:
    const [xc, yc] = xyz.chromaticity();
    assert.assertAlmostEquals(xc, 0.31272);
    assert.assertAlmostEquals(yc, 0.32903);


    // to get the luminous value:
    const l = xyz.luminousValue();
    assert.assertAlmostEquals(l, 100.0);
    // D65 CIE 1931 chromaticity coordinates
    const xyz = new cmt.XYZ(0.31272, 0.32903);
    ```
    */

    #[wasm_bindgen(constructor, variadic)]
    pub fn new_js(x: f64, y:f64, opt : &js_sys::Array) -> Result<XYZ, crate::CmError> {
        use wasm_bindgen::convert::TryFromJsValue;
        use crate::CmError; 
        let (x, y, z, obs) = match opt.length() {
            0 => (x * 100.0/y, 100.0, (1.0 - x - y) * 100.0/y, ObsId::Std1931),
            1 => {
                if opt.get(0).as_f64().is_some() {
                    (x, y, opt.get(0).as_f64().unwrap(), ObsId::Std1931)
                } else {
                    let obs = ObsId::try_from_js_value(opt.get(0))?;
                    (x * 100.0/y, 100.0, (1.0 - x - y) * 100.0/y, obs)
                }
            }
            2 => {
                let z = opt.get(0).as_f64().ok_or(crate::CmError::ErrorString("please provide a z value as number".into()))?;
                let obs = ObsId::try_from_js_value(opt.get(1))?;
                (x, y, z, obs)
            }
            _ => {
                return Err(CmError::ErrorString("Invalid Arguments for XYZ constructor".into()));
            }
        };
        if x<0.0 || y<0.0 || z<0.0 {
                return Err(CmError::ErrorString("XYZ values should be all positive values".into()));
        }
        Ok(XYZ::new(x, y, z, obs))
    }


    /// Get the XYZ tristimulus value as an array.
    #[wasm_bindgen(js_name=values)]
    pub fn values_js(&self)->js_sys::Array {
        let &[x, y, z] = self.data.as_ref();
        js_sys::Array::of3(&x.into(), &y.into(), &z.into())

    }

    /// Get the chromaticity coordinates
    #[wasm_bindgen(js_name=chromaticity)]
    pub fn chromaticity_js(&self)->js_sys::Array {
        let [x, y] = self.chromaticity();
        js_sys::Array::of2(&x.into(), &y.into())
    }

    /// Get the luminous value
    #[wasm_bindgen(js_name=luminousValue)]
    pub fn luminous_value_js(&self)->f64 {
        self.luminous_value()
    }
}