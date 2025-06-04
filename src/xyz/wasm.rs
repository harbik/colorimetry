//! JS-WASM Interface code

use super::XYZ;
use crate::observer::Observer;
use nalgebra::Vector3;
use wasm_bindgen::prelude::wasm_bindgen;

#[wasm_bindgen]
impl XYZ {
    /// Create an XYZ Tristimuls Values object.
    ///
    /// Accepts as arguments
    ///
    /// - x and y chromaticity coordinates only , using the "Cie::Cie1931" observer as default
    /// - x and y chromaticity coordinates, and standard observer ID as 3rd argument
    /// - X, Y, and Z tristimulus values, using the "Cie::Cie1931" observer as default
    /// - X, Y, and Z tristimulus values, and a standard Observer ID as 4th argument
    ///
    /// When only x and y chromaticity coordinates are specified, the luminous
    /// value is set to 100.0 candela per square meter.
    ///
    /// ```javascript, ignore
    /// // Create a new XYZ object using D65 CIE 1931 chromaticity coordinates
    /// const xyz = new cmt.XYZ(0.31272, 0.32903);
    ///
    /// // Get and check the corresponding tristimulus values, with a luminous value
    /// // of 100.0
    /// const [x, y, z] = xyz.values();
    /// assert.assertAlmostEquals(x, 95.047, 5E-3); // D65 wikipedia
    /// assert.assertAlmostEquals(y, 100.0);
    /// assert.assertAlmostEquals(z, 108.883, 5E-3);
    ///
    /// // and get back the orgiinal chromaticity coordinates:
    /// const [xc, yc] = xyz.chromaticity();
    /// assert.assertAlmostEquals(xc, 0.31272);
    /// assert.assertAlmostEquals(yc, 0.32903);
    ///
    /// // to get the luminous value:
    /// const l = xyz.luminousValue();
    /// assert.assertAlmostEquals(l, 100.0);
    /// // D65 CIE 1931 chromaticity coordinates
    /// const xyz = new cmt.XYZ(0.31272, 0.32903);
    /// ```
    #[wasm_bindgen(constructor, variadic)]
    pub fn new_js(x: f64, y: f64, opt: &js_sys::Array) -> Result<XYZ, crate::error::Error> {
        use crate::error::Error;
        use wasm_bindgen::convert::TryFromJsValue;
        let (x, y, z, obs) = match opt.length() {
            0 => (
                x * 100.0 / y,
                100.0,
                (1.0 - x - y) * 100.0 / y,
                Observer::Cie1931,
            ),
            1 => {
                if opt.get(0).as_f64().is_some() {
                    (x, y, opt.get(0).as_f64().unwrap(), Observer::Cie1931)
                } else {
                    let obs = Observer::try_from_js_value(opt.get(0))?;
                    (x * 100.0 / y, 100.0, (1.0 - x - y) * 100.0 / y, obs)
                }
            }
            2 => {
                let z = opt.get(0).as_f64().ok_or(Error::ErrorString(
                    "please provide a z value as number".into(),
                ))?;
                let obs = Observer::try_from_js_value(opt.get(1))?;
                (x, y, z, obs)
            }
            _ => {
                return Err(Error::ErrorString(
                    "Invalid Arguments for XYZ constructor".into(),
                ));
            }
        };
        if x < 0.0 || y < 0.0 || z < 0.0 {
            return Err(Error::ErrorString(
                "XYZ values should be all positive values".into(),
            ));
        }
        Ok(XYZ::from_vecs(Vector3::new(x, y, z), obs))
    }

    /// Get the XYZ tristimulus value as an array.
    #[wasm_bindgen(js_name=values)]
    pub fn values_js(&self) -> js_sys::Array {
        let &[x, y, z] = self.xyz.as_ref();
        js_sys::Array::of3(&x.into(), &y.into(), &z.into())
    }

    /// Get the chromaticity coordinates
    #[wasm_bindgen(js_name=chromaticity)]
    pub fn chromaticity_js(&self) -> js_sys::Array {
        let [x, y] = self.chromaticity().to_array();
        js_sys::Array::of2(&x.into(), &y.into())
    }

    /// Get the luminous value, Y.
    #[wasm_bindgen(js_name=y)]
    pub fn y_js(&self) -> f64 {
        self.y()
    }
}
