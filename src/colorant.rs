use std::{borrow::Cow, ops::{Add, AddAssign, Deref, DerefMut, Mul, MulAssign}};

use colored::Color;
use nalgebra::SVector;

use crate::{gaussian_peak_one, wavelength, wavelengths, CmtError, Filter, Spectrum, NS};


#[derive(Clone, Debug, Default)]
pub struct Colorant(pub(crate) Spectrum);

impl Colorant {

    /// Theoretical spectrum of a perfect grey colorant, consisting of 401
    /// values equal to the value given in the argument, over a range from 380
    /// to 780 nanometer. Mainly used for color mixing calculations.
    pub fn gray(gval: f64) -> Self {
        Self(Spectrum(SVector::<f64,NS>::repeat(gval.clamp(0.0, 1.0))))
    }

    /// Theoretical spectrum of a perfect white colorant, consisting of 401
    /// values over a range from 380 to 780 nanometer. Mainly used for
    /// color mixing calculations.
    pub fn white() -> Self {
        Self::gray(1.0)
    }

    /// Theoretical spectrum of a perfect black color patch, consisting of 401
    /// zero values over a range from 380 to 780 nanometer. Mainly used for
    /// color mixing calculations.
    pub fn black() -> Self {
        Self::gray(0.0)
    }


    /// A Rectangular Band Filter, specified by a central wavelength, and a
    /// width, both in units of meter, or nanometer.
    ///
    /// The filter has a peak value of 1.0
    /// ```rust
    /// # use approx::assert_ulps_eq;
    /// use colorimetry as cmt;
    /// let colorant = cmt::Colorant::top_hat(550.0, 1.0);
    /// let bandfilter: &[f64; cmt::NS] = colorant.as_ref();
    /// assert_ulps_eq!(bandfilter[549-380], 0.0);
    /// assert_ulps_eq!(bandfilter[550-380], 1.0);
    /// assert_ulps_eq!(bandfilter[551-380], 0.0);
    ///
    /// let colorant = cmt::Colorant::top_hat(550.0, 2.0);
    /// let bandfilter: &[f64; cmt::NS] = colorant.as_ref();
    /// assert_ulps_eq!(bandfilter[548-380], 0.0);
    /// assert_ulps_eq!(bandfilter[549-380], 1.0);
    /// assert_ulps_eq!(bandfilter[550-380], 1.0);
    /// assert_ulps_eq!(bandfilter[551-380], 1.0);
    /// assert_ulps_eq!(bandfilter[552-380], 0.0);
    /// 
    /// ```
    pub fn top_hat(center: f64, width: f64) -> Self {
        let [center_m, width_m] = wavelengths([center, width]);
        let left = center_m - width_m/2.0;
        let right = center_m + width_m/2.0;
        let data = SVector::<f64,NS>::from_fn(|i,_j|
            {
              //  let w = (i+380) as f64 * 1E-9;
                let w = wavelength(i+380);
                if w < left - f64::EPSILON || w > right + f64::EPSILON { 0.0}
                else {1.0}
            }
        );
        Self(Spectrum(data))
    }

    /// A Gaussian Filter, specified by a central wavelength, and a
    /// full-width-half-maximum value, both in units of meter, or nanometer.
    ///
    /// The filter has a peak value of 1.0
    pub fn gaussian(center: f64, fwhm: f64) -> Self {
        let [center_m, width_m] = wavelengths([center, crate::sigma_from_fwhm(fwhm)]);
        let data = SVector::<f64,NS>::from_fn(|i,_j|
            gaussian_peak_one((i+380) as f64 * 1E-9, center_m, width_m)
        );
        Self(Spectrum(data))
    }


}

/// Create a Colorant Spectrum from a data slice, with data check.
/// 
/// In this library, the Colorant category represents color filters and color patches, with have spectral values between 0.0 and 1.0.
/// This is the preferred way to crete a colorant, as it does a range check.
impl TryFrom<&[f64]> for Colorant {
    type Error = CmtError;

    fn try_from(data: &[f64]) -> Result<Self, Self::Error> {
        if data.iter().any(|&v|v<0.0 || v>1.0){
            Err(CmtError::OutOfRange { name: "Colorant Spectral Value".into(), low: 0.0, high: 1.0 })
        } else {
            let spectrum = Spectrum::try_from(data)?;
            Ok(Self(spectrum))
        }
    }
}

impl<F> From<F> for Colorant 
    where F: Fn(f64)-> f64 {
    /**
        Colorant from an analytical function, defined over a domain from 0.0 to 1.0, covering the
        wavelength range from 380 to 780 nanometer.

        Values are clamped to a range from 0.0 to 1.0.
        ```rust
        use colorimetry::{Colorant, CIE1931, D65};

        // linear filter from 0.0 to 1.0.
        let tilt: Colorant = (|x:f64|x).into();
        let xy = CIE1931.xyz(&D65, Some(&tilt)).chromaticity();
        approx::assert_abs_diff_eq!(xy.as_ref(), [0.4066, 0.4049].as_ref(), epsilon = 1E-4);

        // parabolic filter
        let parabolic: Colorant = (|x:f64|1.0 - 4.0 * (x - 0.5).powi(2)).into();
        let xy = CIE1931.xyz(&D65, Some(&parabolic)).chromaticity();
        approx::assert_abs_diff_eq!(xy.as_ref(), [0.3466, 0.3862].as_ref(), epsilon = 1E-4);
        ```
    */
    fn from(f: F) -> Self {
        let data = SVector::from_fn(|i,_j|{
            let x = i as f64/(NS-1) as f64;
            f(x).clamp(0.0, 1.0)
        });
        Colorant(Spectrum(data))
    }
}


/// Make colorant data available as a generic [`Filter`] entity, used in particular
/// in the [`crate::Observer`] tristiumulus `xyz`-function.
impl Filter for Colorant {
    fn spectrum(&self) -> Cow<Spectrum> {
        Cow::Borrowed(self)
    }
}

impl Deref for Colorant {
    type Target = Spectrum;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Colorant {

    /// Mutable Access Spectrum methods on references of colorant.
    ///
    /// ```rust
    /// use colorimetry::Colorant;
    /// let mut cth = Colorant::top_hat(500.0, 10.0);
    /// cth.smooth(5.0); // use spectrum's smooth method
    ///
    /// let v = cth[505];
    /// approx::assert_abs_diff_eq!(v, 1.0);
    /// ```
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Add for Colorant {
    type Output = Self;

    /// Adds two colorants
    fn add(self, rhs: Self) -> Self::Output {
        let mut r = self.0 + rhs.0;
        r.clamp(0.0, 1.0);
        Colorant(r)
    }
}

impl Mul<f64> for Colorant {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self(self.0 * rhs)
    }
}

impl Mul<Colorant> for f64 {
    type Output = Colorant;

    fn mul(self, rhs: Colorant) -> Self::Output {
        Colorant(self * rhs.0)
    }
}

impl Mul<Colorant> for Colorant {
    type Output = Self;

    /// Multiplication of two colorants using the `*`-operator.
    /// 
    /// Subtractive Mixing.
    /// ```rust
    /// use colorimetry::Colorant;
    /// let w = Colorant::white();
    /// let b = Colorant::black();
    /// let r: Colorant = w * b;
    /// ```
    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0) // use spectrum multiplication
    }
}

impl Mul<&Colorant> for &Colorant {
    type Output = Colorant;

    /// Multiplication of two colorant references using the `*`-operator.
    /// 
    /// Non-consuming multiplication.
    /// Subtractive Mixing.
    /// ```rust
    /// use colorimetry::Colorant;
    /// let w = Colorant::white();
    /// let b = Colorant::black();
    /// let r: Colorant = &w * &b;
    /// approx::assert_abs_diff_eq(r,b);
    /// ```
    fn mul(self, rhs: &Colorant) -> Self::Output {
        Colorant(&self.0 * &rhs.0) // use spectrum multiplication
    }
}
impl AddAssign<&Self> for Colorant {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += rhs.0 // use spectral multiplication
    }
}