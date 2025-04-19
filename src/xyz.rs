use core::f64;
use std::ops::Add;

use approx::{ulps_eq, AbsDiffEq};
use nalgebra::Vector3;
use crate::{
    geometry::{LineAB, Orientation},
    observer::{self, Observer},
    error::CmtError,
    illuminant::Illuminant,
    rgbspace::RgbSpace,
    spectrum::Spectrum,
    rgb::RGB
};
use wasm_bindgen::prelude::wasm_bindgen; 


const D65A: [f64;3] = [95.04, 100.0, 108.86];
pub const XYZ_D65: XYZ = XYZ::new(&D65A, None, Observer::Std1931);
pub const XYZ_D65WHITE: XYZ = XYZ::new(&D65A, Some(&D65A), Observer::Std1931);


#[wasm_bindgen]
#[derive(Clone, Copy, Debug, PartialEq, Default)]
/// A set of two CIE XYZ Tristimulus values, for a Standard Observer.
/// 
/// One is associated with an illuminant or a reference white value, denoted by the fieldname `xyzn`, and
/// a second, optional value, a set of tristimulus values for a stimulus, or 'a color' value.
/// XYZ values are not often used directly, but form the basis for many colorimetric models, such as
/// CIELAB and CIECAM.
pub struct XYZ {
     pub(crate) observer: Observer,
     pub(crate) xyzn:  Vector3<f64>, // illuminant values
     pub(crate) xyz:  Option<Vector3<f64>>, // stimulus values
}


impl XYZ {

    /// Define a set of [XYZ]-values directly, using an identifier for its
    /// associated observer, such as [Observer::Std1931] or [Observer::Std2015].
    /// 
    pub const fn new(xyzn: &[f64], xyz: Option<&[f64]>, observer: Observer) -> Self {
        let xyzn = Vector3::<f64>::new(xyzn[0], xyzn[1], xyzn[2]);
        let xyz = if let Some(xyz) = xyz {
            Some(Vector3::<f64>::new(xyz[0], xyz[1], xyz[2]))
        } else {
            None
        };
        Self { observer, xyz, xyzn} 
    }

    /// Defines [XYZ] values from [nalgebra::Vector3] values directly.
    pub const fn from_vecs(xyzn: Vector3<f64>, xyz: Option<Vector3<f64>>, observer: Observer) -> XYZ {
        Self { observer, xyz, xyzn} 
    }
    
    /// Create tristimulus values from a chromaticity value, with optional Luminous Value l, and
    /// optional white reference value yw.
    /// 
    /// This produces a illuminant XYZ value, with xyz value set as xyzn.
    pub fn try_from_chromaticity(x: f64, y:f64, l: Option<f64>, observer: Option<Observer>) -> Result<XYZ, CmtError> {
        let l = l.unwrap_or(100.0);
        let observer = observer.unwrap_or_default();
        if (x + y) >= 1.0  {
            Err(CmtError::InvalidChromaticityValues)
        } else {
            let s = l/y;
            let xyz = Vector3::new(x * s, l, (1.0 - x - y)*s);
            Ok(Self::from_vecs(xyz, None, observer))

        }
    }

    pub fn try_from_luv60(u: f64, v: f64, l: Option<f64>, observer: Option<Observer>) ->Result<XYZ, CmtError> {
        let den = 2.0 * u - 8.0 * v + 4.0;
        let x =(3.0 * u)/den;
        let y = (2.0 * v)/den;
        XYZ::try_from_chromaticity(x, y, l, observer)
    }

    /// Try to add two tristimulus values.
    /// Requires sharing the same observer, and no reference white set.
    /// Used for adding illuminants.
    pub fn try_add(&self, other: XYZ) -> Result<XYZ, CmtError> {
        if self.observer == other.observer {
            let data = self.xyzn + other.xyzn;
            match (self.xyz, other.xyz) {
               (None, None) => Ok(XYZ::from_vecs(data, None, self.observer)),
                _ => Err(CmtError::NoReferenceWhiteAllowed)
            }

        } else {
            Err(CmtError::RequireSameObserver)
        }
    }

    /// XYZ Tristimulus values in an an array: [X, Y, Z]
    /// ```
    /// use colorimetry::prelude::*;
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = CIE1931.xyz(&StdIlluminant::D65, None).set_illuminance(100.0);
    /// let [x, y, z] = d65_xyz.into();
    /// // Calculated Spreadsheet Values from CIE Datasets, over a range from 380 to 780nm
    /// assert_ulps_eq!(x, 95.042_267, epsilon = 1E-6);
    /// assert_ulps_eq!(y, 100.0);
    /// assert_ulps_eq!(z, 108.861_036, epsilon = 1E-6);
    /// ```
    pub fn values(&self) -> [f64; 3] {
        (*self).into()
    }

    /// Set the illuminance of an illuminant, either for an illuminant directly,
    /// or for the reference illuminant, in case a color sample XYZ.
    /// ```
    /// use colorimetry::prelude::*;
    /// use approx::assert_ulps_eq;
    /// const D65A: [f64;3] = [95.04, 100.0, 108.86];
    ///
    /// let d65_xyz = CIE1931.xyz(&StdIlluminant::D65, None).set_illuminance(100.0);
    /// assert_ulps_eq!(d65_xyz, XYZ::new(&D65A, None, Observer::Std1931), epsilon = 1E-2);
    /// 
    /// let d65_xyz_sample = CIE1931.xyz(&StdIlluminant::D65, Some(&Colorant::white()));
    /// 
    /// assert_ulps_eq!(d65_xyz_sample, XYZ::new(&D65A, Some(&D65A), Observer::Std1931), epsilon = 1E-2);
    /// ```
    pub fn set_illuminance(mut self, illuminance: f64) -> Self {
        let s = illuminance/self.xyzn.y;
        self.xyzn.iter_mut().for_each(|v|*v = *v * s);
        if let Some(xyz0) = &mut self.xyz {
            // colorant with illuminant
            xyz0.iter_mut().for_each(|v|*v = *v * s)
         };
        self
    }
    
    /// The chromaticity coordinates as an array with an  x and y coordinate
    /// ```
    /// use colorimetry::prelude::*;
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = CIE1931.xyz(&StdIlluminant::D65, None);
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
        let &[_, y, _] = self.xyz.unwrap_or(self.xyzn).as_ref();
        y
    }

    /// CIE 1960 UCS Color Space uv coordinates *Deprecated* by the CIE, but
    /// still used for CCT calculation. Applied to illuminant xyzn values only.
    pub fn uv60(&self) -> [f64; 2] {
        let &[x, y, z] = self.xyz.unwrap_or(self.xyzn).as_ref();
        let den = x + 15.0 * y + 3.0 * z;
        [4.0 * x / den, 6.0 * y / den]
    }
    
    /// The CIE 1964 (U*, V*, W*) color space, also known as CIEUVW, based on
    /// the CIE 1960 UCS.
    /// Still used in Color Rendering Index Calculation.
    /// Uses xyz tristimulus values if present, else uses illuminant's values.
    pub fn uvw64(&self, xyz_ref: XYZ) -> [f64; 3] {
        let yy = self.luminous_value();
        let [ur, vr] = xyz_ref.uv60();
        let [u, v] = self.uv60();
        let ww = 25.0 * yy.powf(1.0 / 3.0) - 17.0;
        let uu = 13.0 * ww * (u - ur);
        let vv = 13.0 * ww * (v - vr);
        [uu, vv, ww]
    }

    /// CIE 1976 CIELUV space, with (u',v') coordinates, calculated for stimulus xyz if present, or else for illuminant.
    pub fn uvprime(&self) -> [f64;2] {
        let &[x, y, z] = self.xyz.unwrap_or(self.xyzn).as_ref();
        let den = x + 15.0 * y + 3.0 * z;
        [4.0 * x / den, 9.0 * y / den]
    }

    //
    pub fn uv_prime_distance(&self, other: &Self) -> f64 {
        let [u1, v1] = self.uvprime();
        let [u2, v2]= other.uvprime();
        (v2-v1).hypot(u2-u1)
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
    pub fn dominant_wavelength(&self, white: XYZ) -> Result<f64, CmtError>{
        let mut sign = 1.0;
        let mut low = self.observer.data().spectral_locus_nm_min();
        let mut high = self.observer.data().spectral_locus_nm_max();
        let mut mid = 540usize;  // 200 fails, as its tail overlaps into the blue region
        if white.observer!=self.observer {
            Err(CmtError::RequireSameObserver)
        } else {
            let [mut x, mut y] = self.chromaticity();
            let [xw, yw] = white.chromaticity();
            // if color point is in the purple rotate it around the white point by 180º, and give wavelength a negative value
            let blue_edge = LineAB::try_new([xw, yw], self.observer.data().spectral_locus_by_nm(low).unwrap().chromaticity()).unwrap();
            let red_edge = LineAB::try_new([xw, yw], self.observer.data().spectral_locus_by_nm(high).unwrap().chromaticity()).unwrap();
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
                let bisect = LineAB::try_new([xw, yw], self.observer.data().spectral_locus_by_nm(mid).unwrap().chromaticity()).unwrap();
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
                let low_ab = LineAB::try_new(white.chromaticity(), self.observer.data().spectral_locus_by_nm(low).unwrap().chromaticity()).unwrap();
                let dlow = low_ab.distance_with_sign(x, y);
                let high_ab = LineAB::try_new(white.chromaticity(), self.observer.data().spectral_locus_by_nm(high).unwrap().chromaticity()).unwrap();
                let dhigh= high_ab.distance_with_sign(x, y);
                if dlow<0.0 || dhigh>0.0 { // not ended up between two lines
                    let s = format!("bisection error in dominant wavelength search:  {dlow} {low} {dhigh} {high}");
                    return Err(CmtError::ErrorString(s));
                }
                let dl = (dlow.abs() * high as f64 + dhigh.abs() * low as f64)/(dlow.abs() + dhigh.abs());
                Ok(sign * dl)
            }
        }
    }


    #[cfg(feature="cct")]
    pub fn cct(self) -> Result<crate::cct::CCT, CmtError> {
        self.try_into()
    }

    /// Convert a set of XYZ tristimulus values to RGB values, using the given RGB space identifier. 
    /// This method requires the luminous value of the reference white, which is typically set to 100.0,
    /// or, less common, 1.0, but any other value can be used as well.
    pub fn rgb(&self, space: Option<RgbSpace>) -> RGB {
        let space = space.unwrap_or_default();
        let xyz = self.xyz.unwrap_or(self.xyzn);
        let ywhite = self.xyzn.y;
        let d = xyz.map(|v| v/ywhite); // normalize to 1.0
        let data = self.observer.data().xyz2rgb(space) * d;
        RGB {
            space,
            observer: self.observer,
            rgb: data,
        }
    }

    pub fn srgb(&self) -> [u8;3] {
        self.rgb(None).into()
    }

}

/*
impl AsRef<[f64;3]> for XYZ {
    fn as_ref(&self) -> &[f64;3] {
        self.xyz.as_ref()
    }
}
 */

/// Get normalized tristimulus values as an array. Normalized to the reference white luminance value yn
/// if present, in case of stimulus values. Else, in case of illuminants only, normalized to
/// a y value of 100.
impl From<XYZ> for [f64;3] {
    fn from(xyz0: XYZ) -> Self {
        let xyz = xyz0.xyz.unwrap_or(xyz0.xyzn);
        let s = 100.0/xyz0.xyzn.y;
        let xyz_array: [f64; 3] = *xyz.as_ref();
        xyz_array.map(|v| v * s)
    }
}

impl AbsDiffEq for XYZ {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let xyz_test = match (self.xyz, other.xyz) {
            (None, None) => true,
            (Some(xyz0), Some(xyz1) ) => 
                xyz0.abs_diff_eq(&xyz1, epsilon),
            _ => false
        };
        self.observer == other.observer &&
        self.xyzn.abs_diff_eq(&other.xyzn, epsilon) &&
        xyz_test
    }
}

impl approx::UlpsEq for XYZ {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        let xyz_test = match (self.xyz, other.xyz) {
            (None, None) => true,
            (Some(xyz0), Some(xyz1) ) => 
                xyz0.ulps_eq(&xyz1, epsilon, max_ulps),
            _ => false
        };
        self.observer == other.observer &&
        self.xyzn.ulps_eq(&other.xyzn, epsilon, max_ulps) &&
        xyz_test
    }
}


impl std::ops::Mul<f64> for XYZ {
    type Output = XYZ;

    /// Multiplication with a right-handed float f64.
    fn mul(mut self, rhs: f64) -> Self::Output {
        if let Some(mut xyz) = self.xyz {
            xyz *= rhs;
        }
        self.xyzn *= rhs;
        self
    }
}

impl std::ops::Mul<XYZ> for f64 {
    type Output = XYZ;

    /// Multiplication of a [`XYZ`]` value on the right of "*" with a float on the left,
    /// resulting in a new [`XYZ`]`.
    fn mul(self, mut rhs: XYZ) -> Self::Output {
        if let Some(mut xyz) = rhs.xyz {
            xyz *= self;
        }
        rhs.xyzn *= self;
        rhs 
    }
}

impl std::ops::Add<XYZ> for XYZ {
    type Output = XYZ;

    /// Add tristimulus values using the "+" operator.
    /// 
    /// Panics if not the same oberver is used.
    /// If either XYZ's has only a xyzn-value, indicating it is an illuminant or a stimulus, the other's xyz-value is requalified
    /// as a direct stimulus, overwriting it's xyz value.
    fn add(mut self, rhs: XYZ) -> Self::Output {
        assert!(self.observer == rhs.observer, "Can not add two XYZ values for different observers");
        (self.xyz, self.xyzn) = match (rhs.xyz, self.xyz) {
            (None, None) =>  (None, self.xyzn + rhs.xyzn),
            (Some(xyz1), Some(xyz2))
                => (Some (xyz1 + xyz2), self.xyzn + rhs.xyzn),
            (Some(xyz), None) => (None, xyz + self.xyzn),
            (None, Some(xyz)) => (None, xyz + rhs.xyzn),
        };
        self
         
    }
}

/*
impl Default for XYZ {
    fn default() -> Self {
        Self { observer: Default::default(), xyz: Default::default(), xyzn: Default::default() }
    }
}
 */

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
    pub fn new_js(x: f64, y:f64, opt : &js_sys::Array) -> Result<XYZ, crate::error::CmtError> {
        use wasm_bindgen::convert::TryFromJsValue;
        use crate::error::CmtError; 
        let (x, y, z, obs) = match opt.length() {
            0 => (x * 100.0/y, 100.0, (1.0 - x - y) * 100.0/y, Observer::Std1931),
            1 => {
                if opt.get(0).as_f64().is_some() {
                    (x, y, opt.get(0).as_f64().unwrap(), Observer::Std1931)
                } else {
                    let obs = Observer::try_from_js_value(opt.get(0))?;
                    (x * 100.0/y, 100.0, (1.0 - x - y) * 100.0/y, obs)
                }
            }
            2 => {
                let z = opt.get(0).as_f64().ok_or(CmtError::ErrorString("please provide a z value as number".into()))?;
                let obs = Observer::try_from_js_value(opt.get(1))?;
                (x, y, z, obs)
            }
            _ => {
                return Err(CmtError::ErrorString("Invalid Arguments for XYZ constructor".into()));
            }
        };
        if x<0.0 || y<0.0 || z<0.0 {
                return Err(CmtError::ErrorString("XYZ values should be all positive values".into()));
        }
        Ok(XYZ::from_vecs(Vector3::new(x, y, z), None, obs))
    }


    /// Get the XYZ tristimulus value as an array.
    /// Values of the stimulus, if present, else the illuminant.
    #[wasm_bindgen(js_name=values)]
    pub fn values_js(&self)->js_sys::Array {
        
        let &[x, y, z] = self.xyz.unwrap_or(self.xyzn).as_ref();
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
    use crate::prelude::*;

    #[test] 
    fn xyz_d65_test(){
        let d65 = CIE1931.xyz_d65();
        let xyz: [f64;3] = d65.into();
        println!("{xyz:?}");
        assert_ulps_eq!(xyz.as_ref(), [95.04, 100.0, 108.867].as_slice(), epsilon = 1E-2);
        let xyz = d65.values();
        assert_ulps_eq!(xyz.as_ref(), [95.04, 100.0, 108.867].as_slice(), epsilon = 1E-2);
    }

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

    #[test]
    fn test_rgb_roundtrip() {
        let rgb_blue = RGB::new(0.0, 0.0, 1.0, Some(Observer::Std1931), Some(RgbSpace::SRGB));
        let xyz_blue = rgb_blue.xyz();
        let xy_blue = xyz_blue.chromaticity();
        assert_ulps_eq!(xy_blue.as_ref(), [0.15, 0.06].as_ref(), epsilon = 1E-5);
        let rgbb = xyz_blue.rgb(None);
        assert_ulps_eq!(rgbb, rgb_blue);
    }

    #[test]
    fn ulps_xyz_test() {
        use approx::assert_ulps_eq;
        use nalgebra::Vector3;
        let xyz0 = XYZ::from_vecs(Vector3::zeros(), None, Observer::Std1931);

        let xyz1 = XYZ::from_vecs(Vector3::new(0.0, 0.0, f64::EPSILON), None,  Observer::Std1931);
        assert_ulps_eq!(xyz0, xyz1, epsilon=1E-5);

        let xyz2 = XYZ::from_vecs(Vector3::new(0.0, 0.0, 2.0 * f64::EPSILON), None,  Observer::Std1931);
        approx::assert_ulps_ne!(xyz0, xyz2);

        // different observer
        #[cfg(feature="supplemental-observers")]
        {
            let xyz3 = XYZ::from_vecs(Vector3::zeros(), None, Observer::Std1964);
            approx::assert_ulps_ne!(xyz0, xyz3);
        }

    }
}
