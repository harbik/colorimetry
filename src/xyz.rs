//! Tristimulus Values
//! ==================
//!
//! The `xyz` module provides core types and operations for working with CIE tristimulus values (X,
//! Y, Z) and chromaticity in Rust.  
//!  
//! This module includes:  
//!  
//! - [`Chromaticity`]: a 2D coordinate type (x, y) representing chromaticity alone.  
//! - [`XYZ`]: a tristimulus struct bundling an (X, Y, Z) vector with its associated [`Observer`].  
//! - Constructors for `XYZ` from raw arrays, chromaticity values, and L\*u\*v coordinates.  
//! - Conversion routines:  
//!   - `XYZ::chromaticity()` to get a `Chromaticity` from an `XYZ` value.  
//!   - `XYZ::rgb()` to convert tristimulus values into out-of-gamut [`WideRgb`] in any `RgbSpace`.  
//!   - `XYZ::dominant_wavelength()` to locate the dominant wavelength relative to a white point.  
//!   - CIE UCS (`uv60`, `uvw64`), CIELUV (`uvprime`, `uv_prime_distance`), and optional CCT.  
//! - Arithmetic and comparison implementations (addition, multiplication, [`AbsDiffEq`], [`UlpsEq`](approx::UlpsEq)) for `XYZ`
//!   that enforce matching observers besided numerical equality.
//! - WASM bindings to expose constructors and accessors for JavaScript.  
//! - A comprehensive test suite (`xyz_test`) covering D65 reference, round-trip conversions,  
//!   gamut checks, and dominant wavelength behavior.  
//!  
//! **Key design points:**  
//!  
//! 1. **Observer-tagged values:** Every `XYZ` value carries an [`Observer`] to ensure  
//!    correct color matching functions and prevent invalid transforms.  
//! 2. **Explicit white reference:** Higher-level models (`CieLab`, `CieCam`, etc.) now  
//!    take the white reference as a separate argument, removing hidden assumptions.  
//! 3. **Out-of-gamut handling:** Conversions to `WideRgb` can produce values outside  
//!    the display gamut; the user can then clamp or compress as needed.  
//! 4. **Physical validity:** While `XYZ` can represent non-physical values, validation  
//!    (e.g., non-negative, within spectral locus) may be added in future refinements.  
//!  

use core::f64;
use std::ops::Add;

use crate::{
    error::Error,
    geometry::{LineAB, Orientation},
    illuminant::Illuminant,
    observer::{self, Observer},
    rgb::RgbSpace,
    rgb::WideRgb,
    spectrum::Spectrum,
};
use approx::{ulps_eq, AbsDiffEq};
use nalgebra::{coordinates::X, ArrayStorage, Vector2, Vector3};
use wasm_bindgen::prelude::wasm_bindgen;

//const D65A: [f64; 3] = [95.04, 100.0, 108.86];
//pub const XYZ_D65_CIE1931: XYZ = XYZ::new(D65A, Observer::Std1931);

/// A chromaticity coordinate with x and y values.
#[wasm_bindgen]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Chromaticity {
    xy: Vector2<f64>,
}

impl Chromaticity {
    /// Returns a new `Chromaticity` object with the given x and y coordinates.
    pub const fn new(x: f64, y: f64) -> Self {
        let xy = Vector2::new(x, y);
        Self { xy }
    }

    /// Returns the x coordinate.
    pub fn x(self) -> f64 {
        self.xy.x
    }

    /// Returns the y coordinate.
    pub fn y(self) -> f64 {
        self.xy.y
    }

    /// Returns the chromaticity coordinates as an array in the format [x, y].
    ///
    /// ```
    /// # use colorimetry::xyz::Chromaticity;
    /// let [x, y] = [0.3127, 0.3290];
    /// let chromaticity = Chromaticity::new(x, y);
    /// assert_eq!(chromaticity.to_array(), [x, y]);
    /// ```
    pub fn to_array(self) -> [f64; 2] {
        *self.xy.as_ref()
    }

    /// Returns the chromaticity coordinate as an `nalgebra` vector.
    pub const fn to_vector(self) -> Vector2<f64> {
        self.xy
    }
}

#[wasm_bindgen]
#[derive(Clone, Copy, Debug, PartialEq, Default)]
/// Represents a color by its tristimulus value XYZ color space.
///
/// The `XYZ` struct represents the tristimulus values (X, Y, Z) and the associated observer.
/// The observer defines the color matching functions used for the conversion.
pub struct XYZ {
    pub(crate) observer: Observer,
    pub(crate) xyz: Vector3<f64>, // tristimulus values
}

impl XYZ {
    /// Creates a new `XYZ` value from raw tristimulus components and an optional observer.
    ///
    /// # Arguments
    /// - `xyz` – A `[f64; 3]` array of tristimulus values `(X, Y, Z)`.  
    ///   These are most often in the 0.0–100.0 range but may lie outside it when derived from
    ///   spectral power distributions.
    /// - `observer` – The color-matching standard observer to use.  
    ///
    /// # Notes
    /// - The `Y` component represents luminous quantity (candela per square meter, lux, or lumen per square meter).
    /// - Allowing values outside 0.0–100.0 lets you handle spectra that produce very bright or non-physical results.
    pub const fn new(xyz: [f64; 3], observer: Observer) -> Self {
        let xyz = Vector3::<f64>::from_array_storage(ArrayStorage([xyz]));
        Self { observer, xyz }
    }

    /// Defines [XYZ] values from [nalgebra::Vector3] values directly.
    pub const fn from_vecs(xyz: Vector3<f64>, observer: Observer) -> XYZ {
        Self { observer, xyz }
    }

    // Creates a new `XYZ` value from chromaticity coordinates, an optional luminous value, and an optional observer.
    ///
    /// # Arguments
    /// - `chromaticity` – A `Chromaticity` struct containing the x and y coordinates.  
    /// - `l` – Optional luminous value (Y). Defaults to `100.0` if `None`, and is used to scale X and Z.  
    /// - `observer` – Optional color-matching standard observer. Defaults to `Observer::Std1931` if `None`.
    ///
    /// # Errors
    /// - Returns `CmtError::InvalidChromaticityValues` if `x + y > 1.0 + ε`, since that cannot represent a valid chromaticity.
    ///
    /// # Notes
    /// - The Y component represents luminous quantity (candela / m², lux, or lumen / m²).  
    /// - Chromaticity-derived XYZ values may lie outside the 0.0–100.0 range when generated from high-intensity or non-physical spectra.
    pub fn from_chromaticity(
        chromaticity: Chromaticity,
        l: Option<f64>,
        observer: Option<Observer>,
    ) -> Result<XYZ, Error> {
        let [x, y] = chromaticity.to_array();
        let l = l.unwrap_or(100.0);
        let observer = observer.unwrap_or_default();

        if (x + y) > 1.0 + f64::EPSILON || x < 0.0 || y < 0.0 {
            Err(Error::InvalidChromaticityValues)
        } else {
            let scale = l / y;
            let xyz = Vector3::new(x * scale, l, (1.0 - x - y) * scale);
            Ok(Self::from_vecs(xyz, observer))
        }
    }

    pub fn from_luv60(
        u: f64,
        v: f64,
        l: Option<f64>,
        observer: Option<Observer>,
    ) -> Result<XYZ, Error> {
        let den = 2.0 * u - 8.0 * v + 4.0;
        let x = (3.0 * u) / den;
        let y = (2.0 * v) / den;
        XYZ::from_chromaticity(Chromaticity::new(x, y), l, observer)
    }

    /// Try to add two tristimulus values.
    /// Requires sharing the same observer, and no reference white set.
    /// Used for adding illuminants.
    pub fn try_add(&self, other: XYZ) -> Result<XYZ, Error> {
        if self.observer == other.observer {
            let data = self.xyz + other.xyz;
            Ok(XYZ::from_vecs(data, self.observer))
        } else {
            Err(Error::RequireSameObserver)
        }
    }

    /// Returns the X value.
    /// ```
    /// use colorimetry::{xyz::XYZ, observer::Observer};
    ///
    /// let xyz = XYZ::new([95.1, 95.0, 27.0], Observer::Std1931);
    /// assert_eq!(xyz.x(), 95.1);
    /// ```
    pub fn x(&self) -> f64 {
        self.values()[0]
    }

    /// Returns the luminous value Y of the tristimulus values.
    ///
    /// This value is associated with different types of photometric quantities -
    /// see Wikipedia's article on [Luminous Intensity](https://en.wikipedia.org/wiki/Luminous_intensity).
    /// As a stimuilus, this value has a unit of candela per square meters,
    /// as an illuminant lux, or lumen per square meter.
    /// ```
    /// use colorimetry::{xyz::XYZ, observer::Observer};
    ///
    /// let xyz = XYZ::new([95.1, 95.0, 27.0], Observer::Std1931);
    /// assert_eq!(xyz.y(), 95.0);
    /// ```
    pub fn y(&self) -> f64 {
        self.values()[1]
    }

    /// Returns the Z value.
    /// ```
    /// use colorimetry::{xyz::XYZ, observer::Observer};
    ///
    /// let xyz = XYZ::new([95.1, 95.0, 27.0], Observer::Std1931);
    /// assert_eq!(xyz.z(), 27.0);
    /// ```
    pub fn z(&self) -> f64 {
        self.values()[2]
    }

    /// Returns the XYZ Tristimulus values in an array on the format [X, Y, Z]
    /// ```
    /// use colorimetry::prelude::*;
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = CIE1931.xyz(&CieIlluminant::D65, None).set_illuminance(100.0);
    /// let [x, y, z] = d65_xyz.values();
    /// // Calculated Spreadsheet Values from CIE Datasets, over a range from 380 to 780nm
    /// assert_ulps_eq!(x, 95.042_267, epsilon = 1E-6);
    /// assert_ulps_eq!(y, 100.0);
    /// assert_ulps_eq!(z, 108.861_036, epsilon = 1E-6);
    /// ```
    pub fn values(&self) -> [f64; 3] {
        *self.xyz.as_ref()
    }

    /// Set the illuminance of an illuminant, either for an illuminant directly,
    /// or for the reference illuminant, in case a color sample XYZ.
    /// ```
    /// use colorimetry::prelude::*;
    /// use approx::assert_ulps_eq;
    /// const D65A: [f64;3] = [95.04, 100.0, 108.86];
    ///
    /// let d65_xyz = CIE1931.xyz(&CieIlluminant::D65, None).set_illuminance(100.0);
    /// assert_ulps_eq!(d65_xyz, XYZ::new(D65A, Observer::Std1931), epsilon = 1E-2);
    ///
    /// let d65_xyz_sample = CIE1931.xyz(&CieIlluminant::D65, Some(&Colorant::white()));
    /// assert_ulps_eq!(d65_xyz_sample, XYZ::new(D65A, Observer::Std1931), epsilon = 1E-2);
    /// ```
    pub fn set_illuminance(mut self, illuminance: f64) -> Self {
        if self.xyz.y > f64::EPSILON && illuminance > f64::EPSILON {
            let s = illuminance / self.xyz.y;
            self.xyz.iter_mut().for_each(|v| *v *= s);
            self
        } else {
            // black override
            XYZ::new([0.0, 0.0, 0.0], self.observer)
        }
    }

    /// Returns the chromaticity coordinates of this `XYZ` value.
    /// ```
    /// use colorimetry::prelude::*;
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = CIE1931.xyz(&CieIlluminant::D65, None);
    /// let chromaticity = d65_xyz.chromaticity();
    /// assert_ulps_eq!(chromaticity.to_array().as_ref(), [0.312_738, 0.329_052].as_slice(), epsilon = 1E-6);
    /// ```
    pub fn chromaticity(&self) -> Chromaticity {
        let [x, y, z] = self.values();
        let s = x + y + z;
        Chromaticity::new(x / s, y / s)
    }

    /// CIE 1960 UCS Color Space uv coordinates *Deprecated* by the CIE, but
    /// still used for CCT calculation. Applied to illuminant xyzn values only.
    pub fn uv60(&self) -> [f64; 2] {
        let &[x, y, z] = self.xyz.as_ref();
        let den = x + 15.0 * y + 3.0 * z;
        [4.0 * x / den, 6.0 * y / den]
    }

    /// The CIE 1964 (U*, V*, W*) color space, also known as CIEUVW, based on
    /// the CIE 1960 UCS.
    /// Still used in Color Rendering Index Calculation.
    /// Uses xyz tristimulus values if present, else uses illuminant's values.
    pub fn uvw64(&self, xyz_ref: XYZ) -> [f64; 3] {
        let yy = self.y();
        let [ur, vr] = xyz_ref.uv60();
        let [u, v] = self.uv60();
        let ww = 25.0 * yy.powf(1.0 / 3.0) - 17.0;
        let uu = 13.0 * ww * (u - ur);
        let vv = 13.0 * ww * (v - vr);
        [uu, vv, ww]
    }

    /// CIE 1976 CIELUV space, with (u',v') coordinates, calculated for stimulus xyz if present, or else for illuminant.
    pub fn uvprime(&self) -> [f64; 2] {
        let &[x, y, z] = self.xyz.as_ref();
        let den = x + 15.0 * y + 3.0 * z;
        [4.0 * x / den, 9.0 * y / den]
    }

    //
    pub fn uv_prime_distance(&self, other: &Self) -> f64 {
        let [u1, v1] = self.uvprime();
        let [u2, v2] = other.uvprime();
        (v2 - v1).hypot(u2 - u1)
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
    pub fn dominant_wavelength(&self, white: XYZ) -> Result<f64, Error> {
        let mut sign = 1.0;
        let wavelength_range = self.observer.data().spectral_locus_wavelength_range();
        let mut low = *wavelength_range.start();
        let mut high = *wavelength_range.end();
        let mut mid = 540usize; // 200 fails, as its tail overlaps into the blue region
        if white.observer != self.observer {
            Err(Error::RequireSameObserver)
        } else {
            let chromaticity = self.chromaticity();
            let [mut x, mut y] = [chromaticity.x(), chromaticity.y()];
            let white_chromaticity = white.chromaticity();
            // if color point is in the purple rotate it around the white point by 180º, and give wavelength a negative value
            let blue_edge = LineAB::new(
                white_chromaticity.to_array(),
                self.observer
                    .data()
                    .xyz_at_wavelength(low)
                    .unwrap()
                    .chromaticity()
                    .to_array(),
            )
            .unwrap();
            let red_edge = LineAB::new(
                white_chromaticity.to_array(),
                self.observer
                    .data()
                    .xyz_at_wavelength(high)
                    .unwrap()
                    .chromaticity()
                    .to_array(),
            )
            .unwrap();
            match (blue_edge.orientation(x, y), red_edge.orientation(x, y)) {
                (Orientation::Colinear, _) => return Ok(380.0),
                (_, Orientation::Colinear) => return Ok(699.0),
                (Orientation::Left, Orientation::Right) => {
                    // mirror point into non-purple region
                    sign = -1.0;
                    x = 2.0 * white_chromaticity.x() - x;
                    y = 2.0 * white_chromaticity.y() - y;
                }
                _ => {} // do nothing
            }
            // start bisectional search
            while high - low > 1 {
                let bisect = LineAB::new(
                    white_chromaticity.to_array(),
                    self.observer
                        .data()
                        .xyz_at_wavelength(mid)
                        .unwrap()
                        .chromaticity()
                        .to_array(),
                )
                .unwrap();
                //   let a = bisect.angle_deg();
                match bisect.orientation(x, y) {
                    Orientation::Left => high = mid,
                    Orientation::Right => low = mid,
                    Orientation::Colinear => {
                        low = mid;
                        high = mid;
                    }
                }
                mid = (low + high) / 2;
            }
            if low == high {
                Ok(sign * low as f64)
            } else {
                let low_ab = LineAB::new(
                    white.chromaticity().to_array(),
                    self.observer
                        .data()
                        .xyz_at_wavelength(low)
                        .unwrap()
                        .chromaticity()
                        .to_array(),
                )
                .unwrap();
                let dlow = low_ab.distance_with_sign(x, y);
                let high_ab = LineAB::new(
                    white.chromaticity().to_array(),
                    self.observer
                        .data()
                        .xyz_at_wavelength(high)
                        .unwrap()
                        .chromaticity()
                        .to_array(),
                )
                .unwrap();
                let dhigh = high_ab.distance_with_sign(x, y);
                if dlow < 0.0 || dhigh > 0.0 {
                    // not ended up between two lines
                    let s = format!("bisection error in dominant wavelength search:  {dlow} {low} {dhigh} {high}");
                    return Err(Error::ErrorString(s));
                }
                let dl = (dlow.abs() * high as f64 + dhigh.abs() * low as f64)
                    / (dlow.abs() + dhigh.abs());
                Ok(sign * dl)
            }
        }
    }

    #[cfg(feature = "cct")]
    pub fn cct(self) -> Result<crate::illuminant::CCT, Error> {
        self.try_into()
    }
    /// Converts a set of **XYZ tristimulus values** to **WideRgb values** using the specified RGB space.
    ///
    /// # Arguments
    ///
    /// - `self`: The XYZ color values to be converted, with a reference white Yn value of 100.
    /// - `rgb_space`: The target RGB space identifier (e.g., `sRGB`, `Adobe RGB`), uses the default
    ///   sRGB space if `None`` is supplied.
    ///
    /// # Returns
    ///
    /// A set of  **WideRgb values**, which can be out of gamut for the specified RGB space.
    /// Use any of the following methods to transform a WideRgb value to a valid RGB value
    ///
    /// - `WideRgb::clamp()`
    /// - `WideRgb::compress()`
    ///
    /// These methods will change the color which will be less saturated, and less bright as the orignial color.
    pub fn rgb(&self, space: Option<RgbSpace>) -> WideRgb {
        let space = space.unwrap_or_default();
        let xyz = self.xyz;
        let d = xyz.map(|v| v / 100.0); // normalize to 1.0
        let data = self.observer.data().xyz2rgb(space) * d;
        WideRgb {
            space,
            observer: self.observer,
            rgb: data,
        }
    }
}

impl From<XYZ> for [f64; 3] {
    /// Converts the tristimulus values to an array on the format [X, Y, Z]
    fn from(xyz: XYZ) -> Self {
        xyz.values()
    }
}

impl AbsDiffEq for XYZ {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.observer == other.observer && self.xyz.abs_diff_eq(&other.xyz, epsilon)
    }
}

impl approx::UlpsEq for XYZ {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.observer == other.observer && self.xyz.ulps_eq(&other.xyz, epsilon, max_ulps)
    }
}

impl std::ops::Mul<f64> for XYZ {
    type Output = XYZ;

    /// Multiplication with a right-handed float f64.
    fn mul(mut self, rhs: f64) -> Self::Output {
        self.xyz *= rhs;
        self
    }
}

impl std::ops::Mul<XYZ> for f64 {
    type Output = XYZ;

    /// Multiplication of a [`XYZ`]` value on the right of "*" with a float on the left,
    /// resulting in a new [`XYZ`]`.
    fn mul(self, mut rhs: XYZ) -> Self::Output {
        rhs.xyz *= self;
        rhs
    }
}

impl std::ops::Add<XYZ> for XYZ {
    type Output = XYZ;

    /// Add tristimulus values using the "+" operator.
    ///
    /// Panics if not the same oberver is used.
    fn add(mut self, rhs: XYZ) -> Self::Output {
        assert!(
            self.observer == rhs.observer,
            "Can not add two XYZ values for different observers"
        );
        self.xyz += rhs.xyz;
        self
    }
}

// JS-WASM Interface code
#[cfg(target_arch = "wasm32")]
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
    pub fn new_js(x: f64, y: f64, opt: &js_sys::Array) -> Result<XYZ, crate::error::Error> {
        use crate::error::Error;
        use wasm_bindgen::convert::TryFromJsValue;
        let (x, y, z, obs) = match opt.length() {
            0 => (
                x * 100.0 / y,
                100.0,
                (1.0 - x - y) * 100.0 / y,
                Observer::Std1931,
            ),
            1 => {
                if opt.get(0).as_f64().is_some() {
                    (x, y, opt.get(0).as_f64().unwrap(), Observer::Std1931)
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
        Ok(XYZ::from_vecs(Vector3::new(x, y, z), None, obs))
    }

    /// Get the XYZ tristimulus value as an array.
    /// Values of the stimulus, if present, else the illuminant.
    #[wasm_bindgen(js_name=values)]
    pub fn values_js(&self) -> js_sys::Array {
        let &[x, y, z] = self.xyz.unwrap_or(self.xyzn).as_ref();
        js_sys::Array::of3(&x.into(), &y.into(), &z.into())
    }

    /// Get the chromaticity coordinates
    #[wasm_bindgen(js_name=chromaticity)]
    pub fn chromaticity_js(&self) -> js_sys::Array {
        let [x, y] = self.chromaticity();
        js_sys::Array::of2(&x.into(), &y.into())
    }

    /// Get the luminous value
    #[wasm_bindgen(js_name=luminousValue)]
    pub fn luminous_value_js(&self) -> f64 {
        self.luminous_value()
    }
}

#[cfg(test)]
mod xyz_test {
    use crate::prelude::*;
    use approx::assert_ulps_eq;

    #[test]
    fn xyz_d65_test() {
        let d65 = CIE1931.xyz_d65();
        let xyz: [f64; 3] = d65.into();
        println!("{xyz:?}");
        assert_ulps_eq!(
            xyz.as_ref(),
            [95.04, 100.0, 108.867].as_slice(),
            epsilon = 1E-2
        );
        let xyz = d65.values();
        assert_ulps_eq!(
            xyz.as_ref(),
            [95.04, 100.0, 108.867].as_slice(),
            epsilon = 1E-2
        );
    }

    #[test]
    fn dominant_wavelength_test() {
        let d65 = CIE1931.xyz_d65().set_illuminance(50.0);

        // 550 nm
        let sl = CIE1931
            .xyz_at_wavelength(550)
            .unwrap()
            .set_illuminance(50.0);
        let t = d65.try_add(sl).unwrap();
        let dl = t.dominant_wavelength(d65).unwrap();
        assert_ulps_eq!(dl, 550.0);

        for wl in 380..=699usize {
            let sl2 = CIE1931.xyz_at_wavelength(wl).unwrap();
            //let [slx, sly] = sl2.chromaticity();
            //println!("sl xy: {slx} {sly}");
            let dl = sl2.dominant_wavelength(d65).unwrap();
            assert_ulps_eq!(dl, wl as f64, epsilon = 1E-10);
        }
    }

    #[test]
    fn dominant_wavelength_purple_test() {
        let d65 = CIE1931.xyz_d65();
        let white_chromaticity = d65.chromaticity();

        // get purple line
        let xyzb = CIE1931.xyz_at_wavelength(380).unwrap();
        let [xb, yb] = xyzb.chromaticity().to_array();
        let xyzr = CIE1931.xyz_at_wavelength(699).unwrap();
        let [xr, yr] = xyzr.chromaticity().to_array();
        let line_t = LineAB::new([xb, yb], [xr, yr]).unwrap();
        for wl in 380..=699usize {
            let sl = CIE1931.xyz_at_wavelength(wl).unwrap();
            let chromaticity = sl.chromaticity();
            let line_u =
                LineAB::new(chromaticity.to_array(), white_chromaticity.to_array()).unwrap();
            let ([xi, yi], t, _) = line_t.intersect(&line_u).unwrap();
            if t > 0.0 && t < 1.0 {
                // see https://en.wikipedia.org/wiki/CIE_1931_color_space#Mixing_colors_specified_with_the_CIE_xy_chromaticity_diagram
                let b = xyzb.set_illuminance(100.0 * (yb * (yr - yi)));
                let r = xyzr.set_illuminance(100.0 * (yr * (yi - yb)));
                let s = b.try_add(r).unwrap();
                let dl = s.dominant_wavelength(d65).unwrap();
                assert_ulps_eq!(dl, -(wl as f64));
            }
        }
    }

    #[test]
    fn test_rgb_roundtrip() {
        let rgb_blue = WideRgb::new(0.0, 0.0, 1.0, Some(Observer::Std1931), Some(RgbSpace::SRGB));
        let xyz_blue = rgb_blue.xyz();
        let xy_blue = xyz_blue.chromaticity().to_array();
        assert_ulps_eq!(xy_blue.as_ref(), [0.15, 0.06].as_ref(), epsilon = 1E-5);
        let rgbb = xyz_blue.rgb(None);
        assert_ulps_eq!(rgbb, rgb_blue);
    }

    #[test]
    fn ulps_xyz_test() {
        use approx::assert_ulps_eq;
        use nalgebra::Vector3;
        let xyz0 = XYZ::from_vecs(Vector3::zeros(), Observer::Std1931);

        let xyz1 = XYZ::from_vecs(Vector3::new(0.0, 0.0, f64::EPSILON), Observer::Std1931);
        assert_ulps_eq!(xyz0, xyz1, epsilon = 1E-5);

        let xyz2 = XYZ::from_vecs(
            Vector3::new(0.0, 0.0, 2.0 * f64::EPSILON),
            Observer::Std1931,
        );
        approx::assert_ulps_ne!(xyz0, xyz2);

        // different observer
        #[cfg(feature = "supplemental-observers")]
        {
            let xyz3 = XYZ::from_vecs(Vector3::zeros(), Observer::Std1964);
            approx::assert_ulps_ne!(xyz0, xyz3);
        }
    }
}
