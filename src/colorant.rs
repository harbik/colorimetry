use std::ops::{Add, AddAssign, Deref, DerefMut, Mul, MulAssign};

use colored::Color;
use nalgebra::SVector;

use crate::{gaussian_peak_one, wavelengths, CmtError, Spectrum, NS};


#[derive(Clone, Debug)]
pub struct Colorant(pub(crate) Spectrum);
impl Colorant {
    /// Create a Colorant Spectrum, with data check.
    /// 
    /// In this library, the Colorant category represents color filters and color patches, with have spectral values between 0.0 and 1.0.
    /// This is the preferred way to crete a colorant, as it does a range check.
    pub fn try_new(data: &[f64]) -> Result<Self, crate::CmtError> {
        if data.iter().any(|&v|v<0.0 || v>1.0){
            Err(CmtError::OutOfRange { name: "Colorant Spectral Value".into(), low: 0.0, high: 1.0 })
        } else {
            let spectrum = Spectrum::try_new(data)?;
            Ok(Self(spectrum))
        }
    }


    /// Theoretical spectrum of a perfect grey colorant, consisting of 401
    /// values equal to the value given in the argument, over a range from 380
    /// to 780 nanometer. Mainly used for color mixing calculations.
    pub fn gray(gval: f64) -> Self {
        Self(Spectrum(SVector::<f64,NS>::repeat(gval.clamp(0.0, 1.0))))
    }

    /// Theoretical spectrum of a perfect white colorant, consisting of 401
    /// 1.0 values over a range from 380 to 780 nanometer. Mainly used for
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
    /// let bandfilter = cmt::Spectrum::band_filter(550.0, 1.0).values();
    /// assert_ulps_eq!(bandfilter[549-380], 0.0);
    /// assert_ulps_eq!(bandfilter[550-380], 1.0);
    /// assert_ulps_eq!(bandfilter[551-380], 0.0);
    ///
    /// let bandfilter = cmt::Spectrum::band_filter(550.0, 2.0).values();
    /// assert_ulps_eq!(bandfilter[548-380], 0.0);
    /// assert_ulps_eq!(bandfilter[549-380], 1.0);
    /// assert_ulps_eq!(bandfilter[550-380], 1.0);
    /// assert_ulps_eq!(bandfilter[551-380], 1.0);
    /// assert_ulps_eq!(bandfilter[552-380], 0.0);
    /// 
    /// ```
    pub fn top_hat(center: f64, width: f64) -> Self {
        let [center_m, width_m] = wavelengths([center, width]);
        let data = SVector::<f64,NS>::from_fn(|i,_j|
            {
                let w = (i+380) as f64 * 1E-9;
                let left = center_m - width_m/2.0;
                let right = center_m + width_m/2.0;
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
    pub fn gaussian(center: f64, width: f64) -> Self {
        let [center_m, width_m] = wavelengths([center, width]);
        let data = SVector::<f64,NS>::from_fn(|i,_j|
            gaussian_peak_one((i+380) as f64 * 1E-9, center_m, width_m)
        );
        Self(Spectrum(data))
    }
}


impl Deref for Colorant {
    type Target = Spectrum;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Colorant {

    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Add for Colorant {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Colorant(self.0 + rhs.0)
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
impl AddAssign<&Self> for Colorant {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += rhs.0
    }
}