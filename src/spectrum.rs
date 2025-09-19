// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

//! Basis for Spectral Data Representations
//! =======================================
//!
//! The field of Colorimetry uses mathematical models to describe the sensations in our mind which we
//! call color.  These models are based the spectral composition of stimuli, essentially rays of light
//! hitting the photosensitive cells in the back of our eyes, and the spectral sensitiviy of these
//! cells. The spectral composition of light, and the objects involved in its processing such as filters
//! and painted patches, is represented by the [Spectrum]-object in this library.
//! The spectral sensitivity of human vision is described by an [`Observer`](crate::observer::Observer).
use core::f64;
use std::{
    borrow::Cow,
    collections::BTreeMap,
    iter::Sum,
    ops::{Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign, RangeInclusive},
};

use approx::AbsDiffEq;

use nalgebra::{DVector, SVector};

use crate::{math::Gaussian, Error};

mod wavelength;
pub use wavelength::{to_wavelength, wavelength, wavelengths};

#[cfg(target_arch = "wasm32")]
mod wasm;

/// The wavelength range of the spectrums supported by this library.
///
/// From 380 to 780 nanometers, inclusive in both ends.
pub const SPECTRUM_WAVELENGTH_RANGE: RangeInclusive<usize> = 380..=780;

/// Number of values in the spectrum. This is 401.
pub const NS: usize = *SPECTRUM_WAVELENGTH_RANGE.end() - *SPECTRUM_WAVELENGTH_RANGE.start() + 1;

/**
This container holds spectral values within a wavelength domain ranging from 380
to 780 nanometers, with an interval size of 1 nanometer and a total of 401
values. It also includes a category tag and an optional 'total' value for the
aggregate value associated with the spectrum.

A `Spectrum` can be constructed from data, but many other construction methods
are available in this library, such as standard illuminants A and D65, Planckian
(Black Body) illuminants, or a `Stimulus` spectrum for a pixel of an sRGB
display.
 */
#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Spectrum(pub(crate) SVector<f64, NS>);

impl Spectrum {
    /// Create a new spectrum from the given data. The data must be the 401 values
    /// from 380 to 780 nm.
    ///
    /// If you have spectral data not exactly matching the range 380-780 nm with 1 nm interval,
    /// you can instead use [`linear_interpolate`](Spectrum::linear_interpolate)
    /// or [`sprague_interpolate`](Spectrum::sprague_interpolate).
    pub const fn new(data: [f64; NS]) -> Self {
        Self(SVector::<f64, NS>::from_array_storage(
            nalgebra::ArrayStorage([data]),
        ))
    }

    /// Constructs a Spectrum from a sparse map of wavelength-value pairs.
    /// The method name reflects that the input is a map from wavelength to value.
    /// The map must contain at least one and at most 401 entries, with keys in the range from 380 to 780
    /// nanometers, inclusive.
    /// # Panics
    /// Panics if the map contains less than 1 or more than 401 entries, or if any key is outside the
    /// range from 380 to 780 nanometers.
    /// # Example
    /// ```rust
    /// use colorimetry::spectrum::Spectrum;
    /// use std::collections::BTreeMap;
    /// let data = ([
    ///     (380, 0.0),
    ///     (400, 0.25),
    ///     (500, 0.5),
    ///     (600, 0.75),
    ///     (700, 1.0),
    /// ]);
    /// let spd = Spectrum::from_wavelength_map(&data);
    /// assert_eq!(spd[380], 0.0);
    /// assert_eq!(spd[400], 0.25);
    /// assert_eq!(spd[500], 0.5);
    /// assert_eq!(spd[600], 0.75);
    /// assert_eq!(spd[700], 1.0);
    /// ```
    pub fn from_wavelength_map(data: &[(usize, f64)]) -> Self {
        if !(1..NS).contains(&data.len()) && data.len() != NS {
            panic!(
                "Need a map with a length between 1 and NS, got a length of {}",
                data.len()
            );
        }
        let mut spd = [0f64; NS];
        for (i, v) in data {
            if *i < *SPECTRUM_WAVELENGTH_RANGE.start() || *i > *SPECTRUM_WAVELENGTH_RANGE.end() {
                panic!(
                    "Wavelength {} out of range, must be between {} and {}",
                    i,
                    SPECTRUM_WAVELENGTH_RANGE.start(),
                    SPECTRUM_WAVELENGTH_RANGE.end()
                );
            }
            spd[*i - SPECTRUM_WAVELENGTH_RANGE.start()] = *v;
        }
        Self(SVector::<f64, NS>::from_array_storage(
            nalgebra::ArrayStorage([spd]),
        ))
    }

    /// This function maps spectral data with irregular intervals or intervals different than 1
    /// nanometer to the standard spectrum as used in this library.
    ///
    /// For domains with a regular interval, the wavelength slice should have a size of two, containing
    /// the minimum and maximum wavelength values, both also in units of meters or nanometers.
    ///
    /// For irregular domains, this function requires a slice of wavelengths and a slice of spectral
    /// data, both of the same size. The wavelengths can be specified in units of meters or nanometers.
    ///
    /// In case of duplicate wavelength values the last data values is used, so it is impossible to
    /// define filters with vertical edges using this method.
    ///
    /// ```rust
    /// // Creates a linear gradient filter, with a zero transmission at 380 nanometer, and full
    /// // transmission at 780 nanometer. This is an example using a uniform wavelength domain as input.
    /// use colorimetry::spectrum::Spectrum;
    /// use approx::assert_ulps_eq;
    /// let data = [0.0, 1.0];
    /// let wl = [380.0, 780.0];
    /// let mut spd = Spectrum::linear_interpolate(&wl, &data).unwrap();
    /// assert_ulps_eq!(spd[380], 0.);
    /// assert_ulps_eq!(spd[380+100], 0.25);
    /// assert_ulps_eq!(spd[380+200], 0.5);
    /// assert_ulps_eq!(spd[380+300], 0.75);
    /// assert_ulps_eq!(spd[380+400], 1.0);
    ///
    /// // Creates a top hat filter, with slanted angles, using an irregular
    /// // wavelength domain.
    /// let data = vec![0.0, 1.0, 1.0, 0.0];
    /// let wl = vec![480.0, 490.0, 570.0, 580.0];
    /// let spd = Spectrum::linear_interpolate(&wl, &data).unwrap();
    /// assert_ulps_eq!(spd[380+0], 0.0);
    /// assert_ulps_eq!(spd[380+100], 0.0);
    /// assert_ulps_eq!(spd[380+110], 1.0);
    /// assert_ulps_eq!(spd[380+190], 1.0);
    /// assert_ulps_eq!(spd[380+200], 0.0);
    /// assert_ulps_eq!(spd[380+300], 0.0);
    /// assert_ulps_eq!(spd[380+400], 0.0);
    /// ```
    pub fn linear_interpolate(wavelengths: &[f64], data: &[f64]) -> Result<Self, Error> {
        let data = match wavelengths.len() {
            2 => linterp(wavelengths.try_into().unwrap(), data)?,
            3.. => linterp_irr(wavelengths, data)?,
            _ => return Err(Error::InterpolateWavelengthError),
        };
        Ok(Self(SVector::<f64, 401>::from_array_storage(
            nalgebra::ArrayStorage([data]),
        )))
    }

    /// Interpolation using Sprague
    ///
    /// This method can only use equally distant data points as input.
    /// See Kerf's paper
    /// [The Interpolation Method of Sprague-Karup](https://www.sciencedirect.com/science/article/pii/0771050X75900273)
    /// for the description of the method.
    /// This implementation uses end-point values for extrapolation, as recommended by CIE15:2004 7.2.2.1.
    pub fn sprague_interpolate(wavelengths: [f64; 2], data: &[f64]) -> Result<Self, Error> {
        let data = sprinterp(wavelengths, data)?;
        Ok(Self(SVector::<f64, 401>::from_array_storage(
            nalgebra::ArrayStorage([data]),
        )))
    }

    pub fn clamp(&mut self, min: f64, max: f64) {
        self.0.iter_mut().for_each(|v| *v = v.clamp(min, max));
    }

    /**
    Smooth a Spectrum by convolution with a Gaussian function
     */
    pub fn smooth(&mut self, mut fwhm: f64) {
        if fwhm < 1E-3 {
            fwhm *= 1E6
        }; // to nanometer
        let gaussian = Gaussian::from_fwhm(0.0, fwhm);
        let sd3 = (6.0 * gaussian.sigma()).floor() as i32;
        let mut kernel = DVector::<f64>::from_iterator(
            (2 * sd3 + 1) as usize,
            (-sd3..=sd3).map(|i| gaussian.peak_one(i as f64)),
        );

        // The smooth operation should not change the energy in a spectrum, so we scale the kernel
        // vector to have a sum of 1.0.
        let sum = kernel.sum();
        kernel.iter_mut().for_each(|v| *v /= sum);

        // use nalgebra's convolve to apply the smooth function, and shift it
        let t = self.0.convolve_full(kernel);
        self.0 = SVector::from_iterator(t.iter().copied().skip(sd3 as usize).take(NS));
    }

    /// Returns the spectral data values, as an array of floats.
    ///
    /// The array contains the 401 data points from 380 to 780 nanometer.
    pub fn values(&self) -> &[f64; NS] {
        &self.0.data.0[0]
    }
}

impl From<[f64; NS]> for Spectrum {
    fn from(data: [f64; NS]) -> Self {
        Self::new(data)
    }
}

impl TryFrom<&[f64]> for Spectrum {
    type Error = Error;

    fn try_from(data: &[f64]) -> Result<Self, Self::Error> {
        if data.len() != NS {
            Err(Error::DataSize401Error)
        } else {
            Ok(Self(SVector::<f64, NS>::from_iterator(
                data.iter().copied(),
            )))
        }
    }
}

impl Default for Spectrum {
    fn default() -> Self {
        Self(SVector::<f64, NS>::zeros())
    }
}

impl AsRef<[f64; 401]> for Spectrum {
    fn as_ref(&self) -> &[f64; 401] {
        self.values()
    }
}

impl AsRef<[f64]> for Spectrum {
    fn as_ref(&self) -> &[f64] {
        &self.0.data.0[0]
    }
}

impl Sum for Spectrum {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut s = Self::default();
        iter.for_each(|si| s += si);
        s
    }
}

impl AbsDiffEq for Spectrum {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

/// Multiplication of two spectra using the `*`-operator, typically for a combinations of an
/// illuminant and a colorant or when combining multiple ColorPatchs or filters.
/// Subtractive Mixing.
impl Mul for Spectrum {
    type Output = Self;

    // multiply two cie spectra
    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0.component_mul(&(rhs.0)))
    }
}

impl Mul<&Spectrum> for &Spectrum {
    type Output = Spectrum;

    /// Multiplication of two references to spectra using the `*`-operator.
    fn mul(self, rhs: &Spectrum) -> Self::Output {
        let s = self.0.component_mul(&(rhs.0));
        Spectrum(s)
    }
}

impl Mul<f64> for Spectrum {
    type Output = Spectrum;

    // spectrum * scalar
    fn mul(self, rhs: f64) -> Self::Output {
        Self(self.0 * rhs)
    }
}

impl Mul<Spectrum> for f64 {
    type Output = Spectrum;

    // scalar * spectrum
    fn mul(self, rhs: Spectrum) -> Self::Output {
        Spectrum(self * rhs.0)
    }
}

impl Mul<&Spectrum> for f64 {
    type Output = Spectrum;

    // scalar * spectrum
    fn mul(self, rhs: &Spectrum) -> Self::Output {
        Spectrum(self * rhs.0)
    }
}

impl Div<&Spectrum> for &Spectrum {
    type Output = Spectrum;

    // multiply two cie spectra
    fn div(self, rhs: &Spectrum) -> Self::Output {
        let s = self.0.component_div(&(rhs.0));
        Spectrum(s)
    }
}

/// Create a Copy On Write instance from a spectrum reference.
impl<'a> From<&'a Spectrum> for Cow<'a, Spectrum> {
    fn from(spectrum: &'a Spectrum) -> Self {
        Cow::Borrowed(spectrum)
    }
}

/// Addition of spectra, typically used for illuminant (multiple sources).
/// Additive mixing
impl Add for Spectrum {
    type Output = Self;

    // add two cie spectra
    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

/// Addition of spectra, typically used for illuminant (multiple sources).
/// Additive mixing
impl Add for &Spectrum {
    type Output = Spectrum;

    // add two cie spectra
    fn add(self, rhs: Self) -> Self::Output {
        let s = self.0 + rhs.0;
        Spectrum(s)
    }
}

/// Addition of spectra, typically used for illuminant (multiple sources).
/// Additive mixing
impl AddAssign for Spectrum {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0
    }
}

/// Addition of spectra, typically used for illuminant (multiple sources).
/// Additive mixing
impl AddAssign<&Spectrum> for Spectrum {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += rhs.0
    }
}

impl MulAssign for Spectrum {
    fn mul_assign(&mut self, rhs: Self) {
        self.0
            .iter_mut()
            .zip(rhs.0.iter())
            .for_each(|(v, w)| *v *= *w);
    }
}

impl MulAssign<f64> for Spectrum {
    fn mul_assign(&mut self, rhs: f64) {
        self.0.iter_mut().for_each(|v| *v *= rhs);
    }
}

/// Read a spectrum value by an integer wavelength value in the range from 380..=780
/// nanometer.
///
/// # Panics
///
/// Panic if the index is less than 380 or greater than 780.
impl Index<usize> for Spectrum {
    type Output = f64;

    fn index(&self, i: usize) -> &Self::Output {
        &self.0[(i - SPECTRUM_WAVELENGTH_RANGE.start(), 0)]
    }
}

/// Mutable Access a spectrum value by an integer wavelength value in the range from 380..=780
/// nanometer.
///
/// # Panics
///
/// Panic if the index is less than 380 or greater than 780.
impl IndexMut<usize> for Spectrum {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.0[(i - SPECTRUM_WAVELENGTH_RANGE.start(), 0)]
    }
}

/// Linear interpolation over a dataset over an equidistant wavelength domain
fn linterp(mut wl: [f64; 2], data: &[f64]) -> Result<[f64; NS], Error> {
    wl.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let [wl, wh] = wavelengths(wl);
    let dlm1 = data.len() - 1; // data length min one

    let mut spd = [0f64; NS];
    spd.iter_mut().enumerate().for_each(|(i, v)| {
        let l = (i + SPECTRUM_WAVELENGTH_RANGE.start()) as f64 * 1E-9; // wavelength in meters
        let t = ((l - wl) / (wh - wl)).clamp(0.0, 1.0); // length parameter
        let tf = t * dlm1 as f64;
        let j = tf.trunc() as usize;
        let f = tf.fract();
        if j >= dlm1 {
            *v = data[dlm1];
        } else {
            *v = data[j] * (1.0 - f) + data[j + 1] * f;
        }
    });
    Ok(spd)
}

/**
Spectrum constructed by linear interpolatino over a dataset with an irregular
wavelength domain.

This algorithm uses a BTreeMap coolection, with wavelengths in picometers as key,
to find a data interval containing the target wavelengths.
 */
fn linterp_irr(wl: &[f64], data: &[f64]) -> Result<[f64; NS], Error> {
    if wl.len() != data.len() {
        Err(Error::InterpolateWavelengthError)
    } else {
        // BTreeMap can not work with floats as keys, using picometer unit
        // (E-12) here as key, so the precision is here three decimals in units
        // of nanometer
        let a = if wl.iter().any(|v| *v > 1E-3) {
            // nanometers
            BTreeMap::from_iter(
                wl.iter()
                    .map(|v| (*v * 1E3) as usize)
                    .zip(data.iter().copied()),
            )
        } else {
            // meters
            BTreeMap::from_iter(
                wl.iter()
                    .map(|v| (*v * 1E12) as usize)
                    .zip(data.iter().copied()),
            )
        };
        let mut spd = [0f64; NS];
        spd.iter_mut().enumerate().for_each(|(i, v)| {
            let k = (i + SPECTRUM_WAVELENGTH_RANGE.start()) * 1000;
            let p = a.range(..k).next_back(); // find values < k
            let n = a.range(k..).next(); // find values >= k
            match (p, n) {
                (Some((&i, &l)), Some((&j, &r))) => {
                    if j == k {
                        *v = r
                    } else {
                        let f = (k - i) as f64 / (j - i) as f64;
                        *v = l * (1.0 - f) + r * f
                    }
                }
                (None, Some((&_j, &r))) => *v = r, // no previous: target wavelength left from lowest value in input dataset, extrapolate
                (Some((&_i, &l)), None) => *v = l, // no next: target wavelength right from highest value in input dataset, extrapolate
                (None, None) => *v = f64::NAN,     // this should never happen
            }
        });
        Ok(spd)
    }
}

/// Sprague interpolation over a dataset over an equidistant wavelength domain
fn sprinterp(mut wl: [f64; 2], data: &[f64]) -> Result<[f64; NS], Error> {
    let imax = data.len() - 1;
    if imax < 6 {
        return Err(Error::ProvideAtLeastNValues(imax));
    };
    let f64_imax = imax as f64;

    // function to deal with extrapolation
    let di = |i: i32| -> f64 {
        if i >= 0 && i <= imax as i32 {
            data[i as usize]
        } else if i < 0 {
            data[0]
        } else {
            data[imax]
        }
    };

    wl.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let [wl, wh] = wavelengths(wl);

    let mut spd = [0f64; NS];
    spd.iter_mut().enumerate().for_each(|(i, v)| {
        let l = (i + SPECTRUM_WAVELENGTH_RANGE.start()) as f64 * 1E-9; // wavelength in meters
        let t = (l - wl) / (wh - wl); // length parameter
        let th = (t * f64_imax).clamp(0.0, f64_imax);
        let h = th.fract();
        let j = th.trunc() as i32;
        *v = sprague(
            h,
            &[di(j - 2), di(j - 1), di(j), di(j + 1), di(j + 2), di(j + 3)],
        );
    });
    Ok(spd)
}

fn sprague(h: f64, v: &[f64]) -> f64 {
    let cf = [
        v[2],
        (v[0] - 8.0 * v[1] + 8.0 * v[3] - v[4]) / 12.0,
        (-v[0] + 16.0 * v[1] - 30.0 * v[2] + 16.0 * v[3] - v[4]) / 24.0,
        (-9.0 * v[0] + 39.0 * v[1] - 70.0 * v[2] + 66.0 * v[3] - 33.0 * v[4] + 7.0 * v[5]) / 24.0,
        (13.0 * v[0] - 64.0 * v[1] + 126.0 * v[2] - 124.0 * v[3] + 61.0 * v[4] - 12.0 * v[5])
            / 24.0,
        (-5.0 * v[0] + 25.0 * v[1] - 50.0 * v[2] + 50.0 * v[3] - 25.0 * v[4] + 5.0 * v[5]) / 24.0,
    ];
    // horner's rule
    cf.into_iter().rev().fold(0.0, |acc, coeff| acc * h + coeff)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::colorant::Colorant;
    use crate::illuminant::{Illuminant, D65};
    use crate::observer::Observer::Cie1931;
    use crate::rgb::Rgb;
    use crate::stimulus::Stimulus;
    use crate::traits::Filter;
    use approx::assert_ulps_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_spectrum_from_rgb() {
        let white: Stimulus = Rgb::new(1.0, 1.0, 1.0, None, None).unwrap().into();
        approx::assert_ulps_eq!(
            Cie1931.xyz_from_spectrum(&white),
            Cie1931.xyz_d65().set_illuminance(100.0),
            epsilon = 5E-4
        );
        let red = Stimulus::from_srgb(255, 0, 0);
        assert_ulps_eq!(
            Cie1931
                .xyz_from_spectrum(&red)
                .chromaticity()
                .to_array()
                .as_ref(),
            &[0.64, 0.33].as_ref(),
            epsilon = 1E-5
        );
    }

    #[test]
    fn test_led() {
        use approx::assert_ulps_eq;
        let ls = Illuminant::led(550.0, 25.0);
        assert_ulps_eq!(ls.irradiance(), 1.0, epsilon = 1E-9);
    }

    #[test]
    fn test_chromaticity() {
        let xyz0 = Cie1931.xyz_from_spectrum(D65.as_ref());
        let [x0, y0] = xyz0.chromaticity().to_array();

        let d65 = D65.clone().set_illuminance(Cie1931, 100.0);
        let xyz = Cie1931.xyz_from_spectrum(d65.as_ref());
        let [x, y] = xyz.chromaticity().to_array();

        assert_ulps_eq!(x0, x);
        assert_ulps_eq!(y0, y);
    }

    #[test]
    fn index_test() {
        let mut s = *Colorant::white().spectrum();

        assert_ulps_eq!(s[500], 1.0);

        // Set a spectral value
        s[500] = 0.5;
        assert_ulps_eq!(s[500], 0.5);

        // Check that neigbhboring values are not changed
        assert_ulps_eq!(s[499], 1.0);
        assert_ulps_eq!(s[501], 1.0);
    }

    #[test]
    #[should_panic]
    fn index_out_of_bounds() {
        let s = *Colorant::white().spectrum();

        let _ = s[300];
    }

    #[test]
    #[should_panic]
    fn index_mut_out_of_bounds() {
        let mut s = *Colorant::white().spectrum();

        s[800] = 0.5;
    }

    #[test]
    fn test_wavelengths() {
        use approx::assert_ulps_eq;

        let mut v1 = [380.0];
        v1 = wavelengths(v1);
        assert_ulps_eq!(v1[0], 380E-9);

        let mut v2 = [380E-9, 780E-9];
        v2 = wavelengths(v2);
        assert_ulps_eq!(v2[0], 380E-9);
        assert_ulps_eq!(v2[1], 780E-9);
    }

    #[test]
    fn ee() {
        let chromaticity = Cie1931
            .xyz_from_spectrum(
                Illuminant::equal_energy()
                    .set_illuminance(Cie1931, 100.0)
                    .as_ref(),
            )
            .chromaticity();
        assert_ulps_eq!(chromaticity.x(), 0.333_3, epsilon = 5E-5);
        assert_ulps_eq!(chromaticity.y(), 0.333_3, epsilon = 5E-5);
    }

    #[test]
    fn d65() {
        let chromaticity = Cie1931
            .xyz_from_spectrum(Illuminant::d65().set_illuminance(Cie1931, 100.0).as_ref())
            .chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(chromaticity.x(), 0.312_72, epsilon = 5E-5);
        assert_ulps_eq!(chromaticity.y(), 0.329_03, epsilon = 5E-5);
    }

    #[test]
    fn d65_test() {
        let [x, y, z] = D65.xyz(None).set_illuminance(100.0).to_array();
        assert_ulps_eq!(x, 95.04, epsilon = 5E-3);
        assert_ulps_eq!(y, 100.0, epsilon = 5E-3);
        assert_ulps_eq!(z, 108.86, epsilon = 5E-3);
    }

    #[test]
    fn d50() {
        let chromaticity = Cie1931
            .xyz_from_spectrum(Illuminant::d50().set_illuminance(Cie1931, 100.0).as_ref())
            .chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(chromaticity.x(), 0.345_67, epsilon = 5E-5);
        assert_ulps_eq!(chromaticity.y(), 0.358_51, epsilon = 5E-5);
    }

    #[test]
    fn add_spectra() {
        use approx::assert_ulps_eq;
        let mut g1 = *Colorant::gray(0.5).spectrum();
        let g2 = *Colorant::gray(0.5).spectrum();
        let g = g1 + g2;
        for i in 380..780 {
            assert_ulps_eq!(g[i], 1.0);
        }

        g1 += &g2;
        for i in 380..780 {
            assert_ulps_eq!(g1[i], 1.0);
        }

        let spectrum = *Colorant::gaussian(550.0, 50.0).spectrum();
        let v = 2.0 * spectrum + -2.0 * spectrum;
        for i in 380..780 {
            assert_ulps_eq!(v[i], 0.0);
        }
    }
    #[test]
    fn mul_spectra_test() {
        use approx::assert_ulps_eq;
        let g = *Colorant::gray(0.5).spectrum();

        let w = 2.0 * g;
        for i in 380..780 {
            assert_ulps_eq!(w[i], 1.0);
        }

        let v = g * 2.0;
        for i in 380..780 {
            assert_ulps_eq!(v[i], 1.0);
        }
    }

    #[test]
    fn test_linterp() {
        use approx::assert_ulps_eq;

        let data = [0.0, 1.0, 0.0];
        let wl = [380.0, 780.0];
        let spd = linterp(wl, &data).unwrap();
        assert_ulps_eq!(spd[0], 0.);
        assert_ulps_eq!(spd[100], 0.5);
        assert_ulps_eq!(spd[200], 1.0);
        assert_ulps_eq!(spd[300], 0.5);
        assert_ulps_eq!(spd[400], 0.0);

        let data = [0.0, 1.0];
        let wl = [380.0, 780.0];
        let spd = linterp(wl, &data).unwrap();
        assert_ulps_eq!(spd[0], 0.);
        assert_ulps_eq!(spd[100], 0.25);
        assert_ulps_eq!(spd[200], 0.5);
        assert_ulps_eq!(spd[300], 0.75);
        assert_ulps_eq!(spd[400], 1.0);

        let data2 = [0.0, 1.0];
        let wl2 = [480.0, 580.0];
        let spd2 = linterp(wl2, &data2).unwrap();
        // print!("{:?}", spd2);
        assert_ulps_eq!(spd2[0], 0.0);
        assert_ulps_eq!(spd2[100], 0.0);
        assert_ulps_eq!(spd2[150], 0.5, epsilon = 1E-10);
        assert_ulps_eq!(spd2[200], 1.0);
        assert_ulps_eq!(spd2[300], 1.0);
        assert_ulps_eq!(spd2[400], 1.0);

        let data3 = [0.0, 1.0];
        let wl3 = [0.0, 1000.0];
        let spd3 = linterp(wl3, &data3).unwrap();
        // print!("{:?}", spd2);
        assert_ulps_eq!(spd3[0], 0.38);
        assert_ulps_eq!(spd3[100], 0.48);
        assert_ulps_eq!(spd3[200], 0.58);
        assert_ulps_eq!(spd3[300], 0.68);
        assert_ulps_eq!(spd3[400], 0.78);
    }

    #[test]
    fn test_smooth() {
        let mut s = *Colorant::default().spectrum();
        s[550] = 1.0;
        s.smooth(5.0);

        let gauss = Gaussian::from_fwhm(550.0, 5.0);
        let sigma = gauss.sigma();

        //  let sigma = sigma_from_fwhm(5.0);
        let w = Colorant::gaussian(550.0, sigma);
        let scale = sigma * (PI * 2.0).sqrt(); // integral of a gaussian
        s.0.iter().zip(w.0 .0.iter()).for_each(|(s, w)| {
            let w = w / scale; // change the reference gaussian colorant to have an integral of 1.0
            approx::assert_abs_diff_eq!(s, &w, epsilon = 1E-8);
        });
    }

    #[test]
    fn sprague_ones() {
        let wl = [380.0, 780.0];
        let data = &[1.0; 81];
        let tinterpolate = sprinterp(wl, data).unwrap();
        tinterpolate
            .iter()
            .for_each(|v| approx::assert_ulps_eq!(*v, 1.0));
    }

    #[test]
    fn sprague_tanh() {
        // Test interpolation of 10nm to 1nm intervals for tanh.
        // This function behaves well w.r.t. constant extrapolation.
        const NF: i32 = 20;
        const NT: i32 = NF * 10;
        let wl = [380.0, 780.0];
        let data: Vec<f64> = (-NF..=NF)
            .map(|i| ((i as f64 / (NF as f64)) * 1.5 * PI).tanh())
            .collect();
        let data_want: Vec<f64> = (-NT..=NT)
            .map(|i| ((i as f64 / (NT as f64)) * 1.5 * PI).tanh())
            .collect();
        let tinterpolate = sprinterp(wl, &data).unwrap();
        tinterpolate
            .iter()
            .zip(data_want.iter())
            .for_each(|(&v, w)| approx::assert_ulps_eq!(v, w, epsilon = 1E-4));
    }

    #[test]
    fn sprague_sin() {
        let wl = [380.0, 780.0];
        let data: Vec<f64> = (0..=80)
            .map(|i| {
                let x = i as f64 / 80.0;
                (x * PI).sin()
            })
            .collect();
        let tinterpolate = sprinterp(wl, &data).unwrap();
        tinterpolate.iter().enumerate().for_each(|(i, &v)| {
            let x = i as f64 / 400.0;
            let y = (x * PI).sin();
            approx::assert_ulps_eq!(y, v, epsilon = 4E-3)
        });
        // non boundary points have very high accuracy
        tinterpolate
            .iter()
            .enumerate()
            .skip(10)
            .take(380)
            .for_each(|(i, &v)| {
                let x = i as f64 / 400.0;
                let y = (x * PI).sin();
                let d = (y - v).abs();
                println!("{i} {y:.4} {v:.4} {d:.6e}");
                approx::assert_ulps_eq!(y, v, epsilon = 5E-10)
            });
    }

    #[test]
    fn test_linterp_irr() {
        use approx::assert_ulps_eq;

        let mut data = vec![0.0, 1.0, 0.0];
        let mut wl = vec![380.0, 480.0, 780.0];
        let mut spd = linterp_irr(&wl, &data).unwrap();
        assert_ulps_eq!(spd[0], 0.);
        assert_ulps_eq!(spd[50], 0.5);
        assert_ulps_eq!(spd[100], 1.0);
        assert_ulps_eq!(spd[250], 0.5);
        assert_ulps_eq!(spd[400], 0.0);

        // top hat with slanted angles
        data = vec![0.0, 1.0, 1.0, 0.0];
        wl = vec![480.0, 490.0, 570.0, 580.0];
        spd = linterp_irr(&wl, &data).unwrap();
        assert_ulps_eq!(spd[0], 0.0);
        assert_ulps_eq!(spd[100], 0.0);
        assert_ulps_eq!(spd[110], 1.0);
        assert_ulps_eq!(spd[190], 1.0);
        assert_ulps_eq!(spd[200], 0.0);
        assert_ulps_eq!(spd[300], 0.0);
        assert_ulps_eq!(spd[400], 0.0);
    }
}
