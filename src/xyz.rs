use core::f64;
use std::ops::{Add, Deref, Mul};

use approx::{assert_ulps_ne, AbsDiffEq, Ulps, UlpsEq};
use nalgebra::{max, Vector3};
use crate::{geometry::{LineAB, Orientation}, obs::ObsId, CmError, CIE1931};
use wasm_bindgen::prelude::wasm_bindgen; 


#[wasm_bindgen]
#[derive(Clone, Copy, Debug, PartialEq)]
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
    
    pub fn from_chromaticity(xy: [f64;2], obs_id: ObsId) -> XYZ {
        let [x, y] = xy;
        Self::new(x, y, 1.0 - x -y, obs_id)
    }

    pub fn try_add(&self, other: XYZ) -> Result<XYZ, CmError> {
        if self.obs_id == other.obs_id {
            let data = self.data + other.data;
            Ok(XYZ {data, obs_id: self.obs_id})
        } else {
            Err(CmError::RequireSameObserver)
        }
    }

    /// XYZ Tristimulus values in a [f64;3] form [X, Y, Z]
    /// ```
    /// use crate::colorimetry::{Spectrum, CIE1931};
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = CIE1931.xyz(&Spectrum::d65_illuminant()).set_illuminance(100.0);
    /// let [x, y, z] = d65_xyz.values();
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
    /// let xy = d65_xyz.chromaticity();
    /// assert_ulps_eq!(xy.as_slice(), [0.312_738, 0.329_052].as_slice(), epsilon = 1E-6);
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

    
    /// The Dominant Wavelength of a color point is the wavelength of spectral
    /// color, obtained from the intersection of a line through a white point
    /// and itself, with the spectral locus.  Points on this line were
    /// historically thought off as having the same hue, but that has been
    /// proven wrong.  This value has limited practical use, but is sometimes
    /// used for color definition of LEDs.
    /// 
    /// The spectral locus, being the boundary of all possible colors in the CIE
    /// 1931 diagram, collapses to one point beyond a wavelength of 699nm. As a
    /// result, the maxium range of dominant wavelengths which can be obtained
    /// is from 380 to 699 nanometer;
    /// 
    pub fn dominant_wavelength(&self, white: XYZ) -> Result<f64, CmError>{
        let mut sign = 1.0;
        let mut low = self.obs_id.observer().spectral_locus_nm_min();
        let mut high = self.obs_id.observer().spectral_locus_nm_max();
        let mut mid = 540usize;  // 200 fails, as its tail overlaps into the blue region
        if white.obs_id!=self.obs_id {
            Err(CmError::RequireSameObserver)
        } else {
            let [mut x, mut y] = self.chromaticity();
            let [xw, yw] = white.chromaticity();
            // if color point is in the purple rotate it around the white point by 180º, and give wavelength a negative value
            let blue_edge = LineAB::try_new([xw, yw], self.obs_id.observer().spectral_locus_by_nm(low).unwrap().chromaticity()).unwrap();
            let red_edge = LineAB::try_new([xw, yw], self.obs_id.observer().spectral_locus_by_nm(high).unwrap().chromaticity()).unwrap();
            match (blue_edge.orientation(x, y), red_edge.orientation(x, y)) {
                (Orientation::Colinear, _) => return Ok(380.0),
                (_, Orientation::Colinear) => return Ok(699.0),
                (Orientation::Left, Orientation::Right) => { // mirror point into non-purple region
                    sign = -1.0;
                    x = 2.0 * xw - x;
                    y = 2.0 * yw - y;
                }
                _ => {} // do nothing
            }
            // start bisectional search
            while high - low > 1 {
                let bisect = LineAB::try_new([xw, yw], self.obs_id.observer().spectral_locus_by_nm(mid).unwrap().chromaticity()).unwrap();
             //   let a = bisect.angle_deg();
                match bisect.orientation(x, y) {
                    Orientation::Left => high = mid,
                    Orientation::Right => low = mid,
                    Orientation::Colinear =>  {
                        low = mid;
                        high = mid;
                    }

                }
                mid = (low + high)/2;
            }
            if low == high {
                Ok(sign * low as f64)
            } else {
                let low_ab = LineAB::try_new(white.chromaticity(), self.obs_id.observer().spectral_locus_by_nm(low).unwrap().chromaticity()).unwrap();
                let dlow = low_ab.distance_with_sign(x, y);
                let high_ab = LineAB::try_new(white.chromaticity(), self.obs_id.observer().spectral_locus_by_nm(high).unwrap().chromaticity()).unwrap();
                let dhigh= high_ab.distance_with_sign(x, y);
                if dlow<0.0 || dhigh>0.0 { // not ended up between two lines
                    let s = format!("bisection error in dominant wavelength search:  {dlow} {low} {dhigh} {high}");
                    return Err(CmError::ErrorString(s));
                }
                let dl = (dlow.abs() * high as f64 + dhigh.abs() * low as f64)/(dlow.abs() + dhigh.abs());
                Ok(sign * dl)
            }
        }
    }


    pub fn cct(&self) -> [f64; 2] {
        todo!()

    }


}

impl AbsDiffEq for XYZ {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.obs_id == other.obs_id && self.data.abs_diff_eq(&other.data, epsilon)
    }
}

impl UlpsEq for XYZ {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.obs_id == other.obs_id && self.data.ulps_eq(&other.data, epsilon, max_ulps)
    }
}

#[test]
fn ulps_xyz_test() {
    use approx::assert_ulps_eq;
    let xyz0 = XYZ::new(0.0, 0.0, 0.0, ObsId::Std1931);

    let xyz1 = XYZ::new(0.0, 0.0, f64::EPSILON, ObsId::Std1931);
    assert_ulps_eq!(xyz0, xyz1);

    let xyz2 = XYZ::new(0.0, 0.0, 2.0*f64::EPSILON, ObsId::Std1931);
    assert_ulps_ne!(xyz0, xyz2);

    let xyz3 = XYZ::new(0.0, 0.0, 0.0, ObsId::Std1976);
    assert_ulps_ne!(xyz0, xyz3);
}

impl Add for XYZ {
    type Output = XYZ;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            data: self.data + rhs.data,
            obs_id: self.obs_id,

        }
    }
}

impl Mul<f64> for XYZ {
    type Output = XYZ;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            data: rhs * self.data,
            obs_id: self.obs_id,

        }
    }
}

impl Mul<XYZ> for f64 {
    type Output = XYZ;

    fn mul(self, rhs: XYZ) -> Self::Output {
        Self::Output {
            data: self * rhs.data,
            obs_id: rhs.obs_id,

        }
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

#[cfg(test)]
mod xyz_test {
    use approx::assert_ulps_eq;
    use crate::{LineAB, CIE1931};

    #[test]
    fn dominant_wavelength_test(){
        let d65 = CIE1931.xyz_d65().set_illuminance(50.0);
        
        // 550 nm
        let sl = CIE1931.spectral_locus_by_nm(550).unwrap().set_illuminance(50.0);
        let t = d65.try_add(sl).unwrap();
        let dl = t.dominant_wavelength(d65).unwrap();
        assert_ulps_eq!(dl, 550.0);

        for wl in 380..=699usize {
            let sl2 = CIE1931.spectral_locus_by_nm(wl).unwrap();
            //let [slx, sly] = sl2.chromaticity();
            //println!("sl xy: {slx} {sly}");
            let dl = sl2.dominant_wavelength(d65).unwrap();
            assert_ulps_eq!(dl, wl as f64, epsilon = 1E-10);

        }


    }

    #[test]
    fn dominant_wavelength_purple_test(){
        let d65 = CIE1931.xyz_d65();
        let [xw, yw] = d65.chromaticity();
        
        // get purple line
        let xyzb = CIE1931.spectral_locus_by_nm(380).unwrap();
        let [xb, yb] = xyzb.chromaticity();
        let xyzr = CIE1931.spectral_locus_by_nm(699).unwrap();
        let [xr, yr] = xyzr.chromaticity();
        let line_t = LineAB::try_new([xb, yb], [xr, yr]).unwrap();
        for wl in 380..=699usize {
            let sl = CIE1931.spectral_locus_by_nm(wl).unwrap();
            let [x, y] = sl.chromaticity();
            let line_u = LineAB::try_new([x, y], [xw, yw]).unwrap();
            let ([xi, yi], t, _) = line_t.intersect(&line_u).unwrap();
            if t>0.0 && t<1.0 {
               // see https://en.wikipedia.org/wiki/CIE_1931_color_space#Mixing_colors_specified_with_the_CIE_xy_chromaticity_diagram
               let b = xyzb.set_illuminance(100.0*(yb * (yr - yi)));
               let r =  xyzr.set_illuminance(100.0*(yr * (yi - yb)));
               let s = b.try_add(r).unwrap();
               let dl = s.dominant_wavelength(d65).unwrap();
               assert_ulps_eq!(dl, -(wl as f64));
            }

        }


    }
}