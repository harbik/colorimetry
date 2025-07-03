//! Tristimulus Values Calculations and Models
//! ==========================================
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

mod chromaticity;
mod dominant;
pub use chromaticity::Chromaticity;

mod rel_xyz;
pub use rel_xyz::RelXYZ;
#[cfg(feature = "gamut-tables")]
pub use rel_xyz::RelXYZGamut;

use core::f64;
use std::fmt::Display;

use crate::{error::Error, observer::Observer, rgb::RgbSpace, rgb::WideRgb};
use approx::AbsDiffEq;
use nalgebra::{ArrayStorage, Vector3};

#[cfg(target_arch = "wasm32")]
mod wasm;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
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
    pub const fn from_vec(xyz: Vector3<f64>, observer: Observer) -> XYZ {
        Self { observer, xyz }
    }

    pub fn is_valid(&self) -> bool {
        if self.xyz.iter().any(|v| *v < 0.0 || !v.is_finite()) {
            return false;
        }
        self.observer
            .spectral_locus()
            .contains(self.chromaticity().to_array())
    }

    // Creates a new `XYZ` value from chromaticity coordinates, an optional luminous value, and an optional observer.
    ///
    /// # Arguments
    /// - `chromaticity` – A `Chromaticity` struct containing the x and y coordinates.  
    /// - `l` – Optional luminous value (Y). Defaults to `100.0` if `None`, and is used to scale X and Z.  
    /// - `observer` – Optional color-matching standard observer. Defaults to `Observer::Cie1931` if `None`.
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
            Ok(Self::from_vec(xyz, observer))
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
            Ok(XYZ::from_vec(data, self.observer))
        } else {
            Err(Error::RequireSameObserver)
        }
    }

    /// Returns the X value.
    /// ```
    /// use colorimetry::{xyz::XYZ, observer::Observer};
    ///
    /// let xyz = XYZ::new([95.1, 95.0, 27.0], Observer::Cie1931);
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
    /// let xyz = XYZ::new([95.1, 95.0, 27.0], Observer::Cie1931);
    /// assert_eq!(xyz.y(), 95.0);
    /// ```
    pub fn y(&self) -> f64 {
        self.values()[1]
    }

    /// Returns the Z value.
    /// ```
    /// use colorimetry::{xyz::XYZ, observer::Observer};
    ///
    /// let xyz = XYZ::new([95.1, 95.0, 27.0], Observer::Cie1931);
    /// assert_eq!(xyz.z(), 27.0);
    /// ```
    pub fn z(&self) -> f64 {
        self.values()[2]
    }

    /// Returns the observer used for this `XYZ` value.
    pub fn observer(&self) -> Observer {
        self.observer
    }

    /// Returns the XYZ Tristimulus values in an array on the format [X, Y, Z]
    /// ```
    /// use colorimetry::prelude::*;
    /// use approx::assert_ulps_eq;
    ///
    /// let d65_xyz = Cie1931.xyz(&CieIlluminant::D65, None).set_illuminance(100.0);
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
    /// let d65_xyz = Cie1931.xyz(&CieIlluminant::D65, None).set_illuminance(100.0);
    /// assert_ulps_eq!(d65_xyz, XYZ::new(D65A, Observer::Cie1931), epsilon = 1E-2);
    ///
    /// let d65_xyz_sample = Cie1931.xyz(&CieIlluminant::D65, Some(&Colorant::white()));
    /// assert_ulps_eq!(d65_xyz_sample, XYZ::new(D65A, Cie1931), epsilon = 1E-2);
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
    /// let d65_xyz = Cie1931.xyz(&CieIlluminant::D65, None);
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

    /// Converts the XYZ tristimulus values to a CCT (Correlated Color Temperature) value.
    ///
    /// # Returns
    /// A `Result` containing the CCT value if the conversion is successful, or an `Error` if it fails.
    ///
    /// # Errors
    /// - `RequiresCIE1931XY`: if the observer is not `CIE1931`. CCT is only defined using the CIE 1931 observer.
    /// - `CCTTemperatureTooHigh`: if `cct` is above 1_000_000 Kelvin.
    /// - `CCTTemperatureTooLow`: if `cct` is below 1000 Kelvin.
    /// - `CCTDuvHighError`: if `duv` is above 0.05.
    /// - `CCTDuvLowError`: if `duv` is below -0.05.
    ///
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
    pub fn rgb(&self, space: RgbSpace) -> WideRgb {
        let xyz = self.xyz;
        let d = xyz.map(|v| v / 100.0); // normalize to 1.0
        let data = self.observer.xyz2rgb_matrix(space) * d;
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

impl std::ops::Sub<XYZ> for XYZ {
    type Output = XYZ;

    /// Subtract tristimulus values using the "-" operator.
    ///
    /// Panics if not the same observer is used.
    fn sub(mut self, rhs: XYZ) -> Self::Output {
        assert!(
            self.observer == rhs.observer,
            "Can not subtract two XYZ values for different observers"
        );
        self.xyz -= rhs.xyz;
        self
    }
}

impl Display for XYZ {
    /// Returns a string representation of the XYZ value in the format "X, Y, Z (Observer)".
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "[{:.3}, {:.3}, {:.3}] ({})",
            self.x(),
            self.y(),
            self.z(),
            self.observer
        )
    }
}

#[cfg(test)]
mod xyz_test {
    use crate::prelude::*;
    use approx::assert_ulps_eq;

    #[test]
    fn xyz_d65_test() {
        let d65 = Cie1931.xyz_d65();
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
    fn test_rgb_roundtrip() {
        use crate::rgb::RgbSpace::SRGB;
        let rgb_blue = WideRgb::new(0.0, 0.0, 1.0, Some(Observer::Cie1931), Some(RgbSpace::SRGB));
        let xyz_blue = rgb_blue.xyz();
        let xy_blue = xyz_blue.chromaticity().to_array();
        assert_ulps_eq!(xy_blue.as_ref(), [0.15, 0.06].as_ref(), epsilon = 1E-5);
        let rgbb = xyz_blue.rgb(SRGB);
        assert_ulps_eq!(rgbb, rgb_blue, epsilon = 1E-6);
    }

    #[test]
    fn ulps_xyz_test() {
        use approx::assert_ulps_eq;
        use nalgebra::Vector3;
        let xyz0 = XYZ::from_vec(Vector3::zeros(), Observer::Cie1931);

        let xyz1 = XYZ::from_vec(Vector3::new(0.0, 0.0, f64::EPSILON), Observer::Cie1931);
        assert_ulps_eq!(xyz0, xyz1, epsilon = 1E-5);

        let xyz2 = XYZ::from_vec(
            Vector3::new(0.0, 0.0, 2.0 * f64::EPSILON),
            Observer::Cie1931,
        );
        approx::assert_ulps_ne!(xyz0, xyz2);

        // different observer
        #[cfg(feature = "supplemental-observers")]
        {
            let xyz3 = XYZ::from_vec(Vector3::zeros(), Observer::Cie1964);
            approx::assert_ulps_ne!(xyz0, xyz3);
        }
    }
}
