/*!
# Spectral Composition

The field of Colorimetry uses mathematical models to describe the sensations in our mind which we
call color.  These models are based the spectral composition of stimuli, essentially rays of light
hitting the photosensitive cells in the back of our eyes, and the spectral sensitiviy of these 
cells. The spectral composition of light, and the objects involved in its processing such as filters
and painted patches, is represented by the [Spectrum]-object in this library.
The spectral sensitivity of human vision is described by an [Observer](crate::Observer).
*/
use core::f64;
use std::{borrow::Cow, collections::BTreeMap, default, error::Error, iter::Sum, ops::{Add, AddAssign, Deref, Div, Index, IndexMut, Mul, MulAssign}};

use num_traits::ToPrimitive;
use url::Url;

use wasm_bindgen::prelude::*;

use nalgebra::{DVector, SVector};

use crate::{
    data::cie_data::{D50, D65, CIE1931},
    observer::ObserverData,
    physics::{gaussian_peak_one, led_ohno, planck, stefan_boltzmann, sigma_from_fwhm, wavelength},
    error::CmtError,
    colorant::Colorant,
    std_illuminants::StdIlluminant,
    physics::C,
    rgb::RGB
};


// Standard Spectrum domain ranging from 380 to 780 nanometer,
// with 401 values.
pub const NS: usize = 401;

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
#[wasm_bindgen]
#[derive(Clone, Copy, Debug)]
pub struct Spectrum(pub(crate) SVector<f64, NS>);

impl Spectrum {


    pub fn mul(mut self, rhs: &Self) -> Self {
        self.0 = self.0.component_mul(&rhs.0) ;
        self
    }

    pub fn mul_f64(mut self, rhs: f64) -> Self {
        self.0 *= rhs;
        self
    }

    /**
    This function maps spectral data with irregular intervals or intervals different than 1
    nanometer to the standard spectrum as used in this library.

    For domains with a regular interval, the wavelength slice should have a size of two, containing
    the minimum and maximum wavelength values, both also in units of meters or nanometers.

    For irregular domains, this function requires a slice of wavelengths and a slice of spectral
    data, both of the same size. The wavelengths can be specified in units of meters or nanometers.

    In case of duplicate wavelength values the last data values is used, so it is impossible to
    define filters with vertical edges using this method.

    ```rust
    // Creates a linear gradient filter, with a zero transmission at 380 nanometer, and full
    // transmission at 780 nanometer. This is an example using a uniform wavelength domain as input.
    use colorimetry as cmt;
    # use approx::assert_ulps_eq;
    let data = [0.0, 1.0];
    let wl = [380.0, 780.0];
    let mut spd = cmt::Spectrum::linear_interpolate(&wl, &data).unwrap().values();
    assert_ulps_eq!(spd[0], 0.);
    assert_ulps_eq!(spd[100], 0.25);
    assert_ulps_eq!(spd[200], 0.5);
    assert_ulps_eq!(spd[300], 0.75);
    assert_ulps_eq!(spd[400], 1.0);

    // Creates a top hat filter, with slanted angles, using an irregular
    // wavelength domain.
    let data = vec![0.0, 1.0, 1.0, 0.0];
    let wl = vec![480.0, 490.0, 570.0, 580.0];
    let spd = cmt::Spectrum::linear_interpolate(&wl, &data).unwrap().values();
    assert_ulps_eq!(spd[0], 0.0);
    assert_ulps_eq!(spd[100], 0.0);
    assert_ulps_eq!(spd[110], 1.0);
    assert_ulps_eq!(spd[190], 1.0);
    assert_ulps_eq!(spd[200], 0.0);
    assert_ulps_eq!(spd[300], 0.0);
    assert_ulps_eq!(spd[400], 0.0);
    ```
    */
    pub fn linear_interpolate(wavelengths: &[f64], data: &[f64]) ->Result<Self, CmtError> {
        let data = match wavelengths.len() {
           2 =>  linterp(wavelengths.try_into().unwrap(), data)?,
           3.. => linterp_irr(wavelengths, data)?,
           _ => return Err(CmtError::InterpolateWavelengthError)
        };
        Ok(Self(SVector::<f64, 401>::from_array_storage(nalgebra::ArrayStorage([data]))))
    }

    
    /// Interpolation using Sprague
    /// 
    /// This method can only use equally distant data points as input.
    /// See Kerf's paper
    /// [The Interpolation Method of Sprague-Karup](https://www.sciencedirect.com/science/article/pii/0771050X75900273)
    /// for the description of the method.
    /// This implementation uses end-point values for extrapolation, as recommended by CIE15:2004 7.2.2.1.
    
    pub fn sprague_interpolate(wavelengths: [f64;2], data: &[f64]) ->Result<Self, CmtError> {
        let data = sprinterp(wavelengths.try_into().unwrap(), data)?;
        Ok(Self(SVector::<f64, 401>::from_array_storage(nalgebra::ArrayStorage([data]))))
    }

    pub fn clamp(&mut self, min: f64, max: f64) {
        self.0.iter_mut().for_each(|v|*v = v.clamp(min, max));
    }


    /**
    Smooth a Spectrum by convolution with a Gaussian function
     */
    pub fn smooth(&mut self, mut fwhm: f64) {
        if fwhm < 1E-3 { fwhm *= 1E6 }; // to nanometer
        let sigma = sigma_from_fwhm(fwhm);
        let sd3 = (6.0 * sigma).floor() as i32;
        let mut kernel =  
            DVector::<f64>
                ::from_iterator(
                    (2*sd3+1) as usize,
                    (-sd3..=sd3)
                        .into_iter()
                        .map(|i| gaussian_peak_one(i as f64, 0.0, sigma)
                    ));

        // The smooth operation should not change the energy in a spectrum, so we scale the kernel
        // vector to have a sum of 1.0.
        let sum = kernel.sum();
        kernel.iter_mut().for_each(|v|*v /= sum);

        // use nalgebra's convolve to apply the smooth function, and shift it
        let t = self.0.convolve_full(kernel);
        self.0 = SVector::from_iterator(t.iter().copied().skip(sd3 as usize).take(NS));
    }



    /// Downloads a spectrum
    pub async fn fetch(loc: &str) -> Result<Self, Box<dyn Error>> {
        let _url = Url::parse(loc)?;
        todo!()


    }


}

impl TryFrom<&[f64]> for Spectrum {
    type Error = CmtError;

    fn try_from(data: &[f64]) -> Result<Self, Self::Error> {
        if data.len()!=NS {
            Err(CmtError::DataSize401Error)
        } else {
            Ok(Self(SVector::<f64, NS>::from_iterator(data.into_iter().copied())))
        }
    }
}

impl Default for Spectrum {
    fn default() -> Self {
        Self (SVector::<f64, NS>::zeros()) 
    }
}

impl AsRef<[f64;401]> for Spectrum {
    fn as_ref(&self) -> &[f64;401] {
        &self.0.data.0[0]
    }
}

impl AsRef<[f64]> for Spectrum {
    fn as_ref(&self) -> &[f64] {
        &self.0.data.0[0]
    }
}

impl Sum for Spectrum {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut s = Self::default() ;
        iter.for_each(|si|s += si);
        s
    }
}

// JS-WASM Interface code
#[cfg(target_arch="wasm32")]
#[wasm_bindgen]
impl Spectrum {

    /// Creates a new Spectrum object, using as input a `Category`, a
    /// Float64Array with exactly 401 datapoints, and an optional third
    /// parameter called total, representing the total irradiance, transmission,
    /// or reflectivity of the values, depending on the category of the
    /// spectrum. The spectral values should be associated with a wavelength
    /// domain from 380 to 480 nanometer, with an interval size of 1 nanometer.
    ///
    /// If the Spectral data you have uses another wavelength domain and/or a different
    /// wavelength interval, use the linear or sprague interpolate constructors,
    /// which takes a wavelength domain and spectral data as arguments.
    #[wasm_bindgen(constructor)]
    pub fn new_js(data: &[f64]) -> Result<Spectrum, wasm_bindgen::JsError> {
        Ok(Spectrum::try_from(data)?)
    }

    /// Returns the spectral data values, as a Float64Array containing 401 data
    /// points, over a wavelength domain from 380 t0 780 nanometer, with a
    /// stepsize of 1 nanometer.
    #[wasm_bindgen(js_name=Values)]
    pub fn values_js(&self) -> Box<[f64]> {
        let values: &[f64] = self.as_ref();
        values.into()
    }

    /**
    This function maps spectral data with irregular intervals or intervals
    different than 1 nanometer to the standard spectrum as used in this
    library.

    For domains with a regular interval, the wavelength slice should have a size
    of two, containing the minimum and maximum wavelength values, both also in
    units of meters or nanometers.

    For irregular domains, this function requires a slice of wavelengths and
    a slice of spectral data, both of the same size. The wavelengths can be
    specified in units of meters or nanometers.

    In case of duplicate wavelength values the last data values is used, so it
    is impossible to define filters with vertical edges using this method.

    ```ts, ignore
    // Creates a linear gradient filter, with a zero transmission at 380
    // nanometer, and full transmission at 780 nanometer. This is an example
    // using a uniform wavelength domain as input.
    use colorimetry as cmt;
    # use approx::assert_ulps_eq;
    let data = [0.0, 1.0];
    let wl = [380.0, 780.0];
    let mut spd = cmt::Spectrum::linear_interpolate(cmt::Category::Colorant, &wl, &data, None).unwrap().values();
    assert_ulps_eq!(spd[0], 0.);
    assert_ulps_eq!(spd[100], 0.25);
    assert_ulps_eq!(spd[200], 0.5);
    assert_ulps_eq!(spd[300], 0.75);
    assert_ulps_eq!(spd[400], 1.0);

    // Creates a top hat filter, with slanted angles, using an irregular
    // wavelength domain.
    let data = vec![0.0, 1.0, 1.0, 0.0];
    let wl = vec![480.0, 490.0, 570.0, 580.0];
    let spd = cmt::Spectrum::linear_interpolate(cmt::Category::Colorant, &wl, &data, None).unwrap().values();
    assert_ulps_eq!(spd[0], 0.0);
    assert_ulps_eq!(spd[100], 0.0);
    assert_ulps_eq!(spd[110], 1.0);
    assert_ulps_eq!(spd[190], 1.0);
    assert_ulps_eq!(spd[200], 0.0);
    assert_ulps_eq!(spd[300], 0.0);
    assert_ulps_eq!(spd[400], 0.0);
    ```
    */
    #[wasm_bindgen(js_name=linearInterpolate)]
    pub fn linear_interpolate_js(wavelengths: &[f64], data: &[f64], total_js: &JsValue) -> Result<Spectrum, CmtError> {
        Self::linear_interpolate(wavelengths, data)

    }

    /// Calculates the Color Rendering Index values for illuminant spectrum.
    /// 
    /// To use this function, first use `await CRI.init()`, which downloads the
    /// Test Color Samples required for the calculation.  These are downloaded
    /// seperately to limit the size of the main web assembly library.
    #[cfg(feature="cri")]
    #[wasm_bindgen(js_name=cri)]
    pub fn cri_js(&self) -> Result<crate::cri::CRI, CmtError> {
        todo!()
    }

}


// Multiplication of two spectra using the `*`-operator, typically for a combinations of an illuminant and a colorant
// or when combining multiple ColorPatchs or filters. Subtractive Mixing.
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
impl <'a> From<&'a Spectrum> for Cow<'a, Spectrum> {
    
    fn from(spectrum: &'a Spectrum) -> Self {
        Cow::Borrowed(spectrum)
    }
}


// Addition of spectra, typically used for illuminant (multiple sources).
// Additive mixing
impl Add for Spectrum {
    type Output = Self;

    // add two cie spectra
    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

// Addition of spectra, typically used for illuminant (multiple sources).
// Additive mixing
impl Add for &Spectrum {
    type Output = Spectrum;

    // add two cie spectra
    fn add(self, rhs: Self) -> Self::Output {
        let s = self.0 + rhs.0;
        Spectrum(s)
    }
}

// Addition of spectra, typically used for illuminant (multiple sources).
// Additive mixing
impl AddAssign for Spectrum {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0
    }
}

// Addition of spectra, typically used for illuminant (multiple sources).
// Additive mixing
impl AddAssign<&Spectrum> for Spectrum {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += rhs.0
    }
}

impl MulAssign for Spectrum {
    fn mul_assign(&mut self, rhs: Self) {
        self.0.iter_mut().zip(rhs.0.iter()).for_each(|(v,w)| *v *= *w);

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
/// Invalid indices will result in a f64::NAN value.
impl Index<usize> for Spectrum {
    type Output = f64;

    fn index(&self, i: usize) -> &Self::Output {
        if i<380 || i>780 {
            &f64::NAN
        } else {
            &self.0[(i-380,0)]

        }
    }
}

/// Mutable Access a spectrum value by an integer wavelength value in the range from 380..=780
/// nanometer.
/// 
/// Out of range indices will set the edge values.
impl IndexMut<usize> for Spectrum {

    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        match i {
            380..=780 => &mut self.0[(i-380,0)],
            ..380 => &mut self.0[(0,0)],
            781.. => &mut self.0[(400,0)],
        }
    }
}

/// Convenience function for specifying wavelengths in nanometers or meters.
///
/// This accepts integer and float values.
/// Wwavelength values larger than 1E-3 are assumed to have the unit nanometer
/// and are converted to a unit of meters.
/// All integer values are nanometaer values.
pub fn wavelengths<T: ToPrimitive, const N: usize>(v:[T; N]) -> [f64;N] {
    v.map(|x|wavelength(x))
}

/// Linear interpolatino over a dataset over an equidistant wavelength domain
fn linterp(mut wl: [f64;2], data: &[f64]) -> Result<[f64;NS], CmtError> {
    wl.sort_by(|a, b| a.partial_cmp(b).unwrap()); 
    let [wl, wh] = wavelengths(wl);
    let dlm1 = data.len()-1; // data length min one
    
    let mut spd = [0f64; NS];
    spd.iter_mut().enumerate().for_each(|(i,v)|{
        let l = (i + 380) as f64 * 1E-9; // wavelength in meters
        let t = ((l-wl)/(wh - wl)).clamp(0.0, 1.0); // length parameter
        let tf = (t * dlm1 as f64) as f64;
        let j = tf.trunc() as usize;
        let f = tf.fract();
        if j >= dlm1 {
            *v = data[dlm1];
        } else {
            *v = data[j] * (1.0 - f) + data[j+1] * f;
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
fn linterp_irr(wl: &[f64], data: &[f64]) -> Result<[f64;NS], CmtError> {
    if wl.len()!=data.len() {
        Err(CmtError::InterpolateWavelengthError)
    } else {
        // BTreeMap can not work with floats as keys, using picometer unit
        // (E-12) here as key, so the precision is here three decimals in units
        // of nanometer
        let a = 
            if wl.iter().any(|v|*v>1E-3) { // nanometers 
                BTreeMap::from_iter(wl.iter().map(|v|(*v*1E3) as usize).zip(data.iter().copied()))
            } else { // meters
                BTreeMap::from_iter(wl.iter().map(|v|(*v*1E12) as usize).zip(data.iter().copied()))
            };
        let mut spd = [0f64; NS];
        spd.iter_mut().enumerate().for_each(|(i,v)|{
            let k = (i + 380) * 1000;
            let p = a.range(..k).next_back(); // find values < k
            let n = a.range(k..).next(); // find values >= k
            match (p,n) {
                (Some((&i, &l)), Some((&j, &r))) => {
                    if j == k { *v = r}
                    else {
                        let f = (k - i) as f64 / (j-i) as f64;
                        *v = l * (1.0 - f) + r * f
                    }
                }
                (None, Some((&_j, &r))) => *v = r, // no previous: target wavelength left from lowest value in input dataset, extrapolate 
                (Some((&_i, &l)), None) => *v = l, // no next: target wavelength right from highest value in input dataset, extrapolate
                (None, None) => *v = f64::NAN // this should never happen
            }
        });
        Ok(spd)
        
    }
}

/// Sprague interpolation over a dataset over an equidistant wavelength domain
fn sprinterp(mut wl: [f64;2], data: &[f64]) -> Result<[f64;NS], CmtError> {
    let imax = data.len()-1;
    if imax <6 { return Err(CmtError::ProvideAtLeastNValues(imax))};
    let f64_imax = imax as f64;

    // function to deal with extrapolation
    let di = |i:i32| -> f64 {
        if i>=0 && i<=imax as i32 { data[i as usize] }
        else if i<0 { data[0] } 
        else { data[imax] }
    };

    wl.sort_by(|a, b| a.partial_cmp(b).unwrap()); 
    let [wl, wh] = wavelengths(wl);

    let mut spd = [0f64; NS];
    spd.iter_mut().enumerate().for_each(|(i,v)|{
        let l = (i + 380) as f64 * 1E-9; // wavelength in meters
        let t = (l-wl)/(wh - wl); // length parameter
        let th = (t * f64_imax).clamp(0.0, f64_imax);
        let h = th.fract();
        let j = th.trunc() as i32;
        *v = sprague(h, &[di(j-2), di(j-1), di(j), di(j+1), di(j+2), di(j+3)]);
    });
    Ok(spd)
}

fn sprague(h: f64, v: &[f64]) -> f64 {
    let cf = [
        v[2],
        (v[0] - 8.0 * v[1] + 8.0 * v[3] - v[4]) / 12.0,
        (-v[0] + 16.0 * v[1] - 30.0 * v[2] + 16.0 * v[3] - v[4]) / 24.0,
        (-9.0 * v[0] + 39.0 * v[1] - 70.0 * v[2] + 66.0 * v[3] - 33.0 * v[4] + 7.0 * v[5]) / 24.0,
        (13.0 * v[0] - 64.0 * v[1] + 126.0 * v[2] - 124.0 * v[3] + 61.0 * v[4] - 12.0 * v[5]) / 24.0,
        (-5.0 * v[0] + 25.0 * v[1] - 50.0 * v[2] + 50.0 * v[3] - 25.0 * v[4] + 5.0 * v[5]) / 24.0,
    ];
    // horner's rule
    cf.into_iter().rev().fold(0.0, |acc, coeff| acc * h + coeff)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use approx::assert_ulps_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_spectrum_from_rgb(){
        let white: Stimulus = RGB::new(1.0, 1.0, 1.0, None, None).into();
        approx::assert_ulps_eq!(CIE1931.xyz_from_spectrum(&white, None), CIE1931.xyz_d65().set_illuminance(100.0), epsilon = 1E-6);
        let red = Stimulus::srgb(255, 0, 0);
        assert_ulps_eq!(CIE1931.xyz_from_spectrum(&red, None).chromaticity().as_ref(), &[0.64, 0.33].as_ref(), epsilon = 1E-5);
    }

    #[test]
    fn test_led(){
        use approx::assert_ulps_eq;
        let ls = Illuminant::led(550.0, 25.0);
        assert_ulps_eq!(ls.irradiance(), 1.0, epsilon = 1E-9);
    }

    #[test]
    fn test_chromaticity(){
        let xyz0 = CIE1931.xyz_from_spectrum(&D65, None);
        let [x0, y0] = xyz0.chromaticity();

        let illuminance = D65.illuminance(&CIE1931);
        let d65 = D65.clone().set_illuminance(&CIE1931, 100.0);
        let xyz = CIE1931.xyz_from_spectrum(&d65, None);
        let [x, y] = xyz.chromaticity();
    
        assert_ulps_eq!(x0, x);
        assert_ulps_eq!(y0, y);
    }

    #[test]
    fn index_test(){
        let mut s = Colorant::white();

        // Set a spectral value
        s[500] = 0.5;
        assert_ulps_eq!(s[500], 0.5);

        // A read from an index out of the 380..=780 range with given f64::NAN.
        s[300] = 0.5;
        assert!(s[300].is_nan());
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
        let [x, y ] = CIE1931.xyz_from_spectrum(
            &Illuminant::equal_energy().set_illuminance(&CIE1931, 100.0), None).chromaticity();
        assert_ulps_eq!(x, 0.333_3, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.333_3, epsilon = 5E-5);
    }

    #[test]
    fn d65() {
        let [x, y ] = CIE1931.xyz_from_spectrum(
            &&Illuminant::d65().set_illuminance(&CIE1931, 100.0), None).chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(x, 0.312_72, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.329_03, epsilon = 5E-5);
    }

    #[test]
    fn d50() {
        let [x, y ] = CIE1931.xyz_from_spectrum(&Illuminant::d50().set_illuminance(&CIE1931, 100.0), None).chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(x, 0.345_67, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.358_51, epsilon = 5E-5);
    }

    #[cfg_attr(test, cfg(feature="cie-illuminants"))]
    fn a() {
        let a: Illuminant = StdIlluminant::A.into();
        let [x, y ] = CIE1931.xyz_from_spectrum(&a, None).chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(x, 0.447_58, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.407_45, epsilon = 5E-5);
    }

    #[test]
    fn add_spectra(){
        use approx::assert_ulps_eq;
        let mut g1 = Colorant::gray(0.5);
        let g2 = Colorant::gray(0.5);
        let g = g1.clone() + g2.clone();
        for i in 380..780 {
            assert_ulps_eq!(g[i], 1.0);
        }

        g1 += &g2;
        for i in 380..780 {
            assert_ulps_eq!(g1[i], 1.0);
        }

        let v = 2.0 * *Colorant::gaussian(550.0, 50.0) + -2.0 * *Colorant::gaussian(550.0, 50.0);
        for i in 380..780 {
            assert_ulps_eq!(v[i], 0.0);
        }

    }
    #[test]
    fn mul_spectra_test(){
        use approx::assert_ulps_eq;
        let g = Colorant::gray(0.5);
    
        let w = 2.0 * g.clone();
        for i in 380..780 {
            assert_ulps_eq!(w[i], 1.0);
        }
    
        let v = g * 2.0;
        for i in 380..780 {
            assert_ulps_eq!(v[i], 1.0);
        }
    
    }

    #[test]
    fn test_linterp(){
        use approx::assert_ulps_eq;

        let data = [0.0, 1.0,  0.0];
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
        let mut s = Colorant::default();
        s[550] = 1.0;
        s.smooth(5.0);

        let w = Colorant::gaussian(550.0, 5.0);

        s.0.0.iter().zip(w.0.0.iter()).enumerate().for_each(|(i, (s,w))|{
            let j = i + 380;
           // println!("{j} {s:.6} {w:.6}");
            approx::assert_abs_diff_eq!(s,w, epsilon=1E-8);
        });

    }

    #[test]
    fn sprague_ones() {
        let wl = [380.0, 780.0];
        let data = &[1.0; 81];
        let tinterpolate = sprinterp(wl, data).unwrap();
        tinterpolate.iter().for_each(|v|approx::assert_ulps_eq!(*v, 1.0));
    }

    #[test]
    fn sprague_tanh() {
        // Test interpolation of 10nm to 1nm intervals for tanh.
        // This function behaves well w.r.t. constant extrapolation.
        const NF: i32 = 20;
        const NT: i32 = NF * 10;
        let wl = [380.0, 780.0];
        let data: Vec<f64> = (-NF..=NF).into_iter().map(|i|((i as f64/(NF as f64))*1.5*PI).tanh()).collect();
        let data_want: Vec<f64> = (-NT..=NT).into_iter().map(|i|((i as f64/(NT as f64))*1.5*PI).tanh()).collect();
        let tinterpolate = sprinterp(wl, &data).unwrap();
        tinterpolate.iter().zip(data_want.iter()).for_each(|(&v, w)|approx::assert_ulps_eq!(v, w, epsilon=1E-4));
    }

    #[test]
    fn sprague_sin() {
        let wl = [380.0, 780.0];
        let data: Vec<f64> = (0..=80).into_iter().map(|i|{
            let x = i as f64/80.0;
            (x*PI).sin()
        }).collect();
        let tinterpolate = sprinterp(wl, &data).unwrap();
        tinterpolate.iter().enumerate().for_each(|(i, &v)|{
            let x = i as f64/400.0;
            let y = (x*PI).sin();
            let d = (y-v).abs();
          //  println!("{i} {y:.4} {v:.4} {d:.6}");
            approx::assert_ulps_eq!(y,v, epsilon=4E-3)
        });
        // non boundary points have very high accuracy
        tinterpolate.iter().enumerate().skip(10).take(380).for_each(|(i, &v)|{
            let x = i as f64/400.0;
            let y = (x*PI).sin();
            let d = (y-v).abs();
            println!("{i} {y:.4} {v:.4} {d:.6e}");
            approx::assert_ulps_eq!(y,v, epsilon=5E-10)
        });
    }

    #[test]
    fn test_linterp_irr(){
        use approx::assert_ulps_eq;

        let mut data = vec![0.0, 1.0, 0.0];
        let mut wl = vec![380.0, 480.0,  780.0];
        let mut spd = linterp_irr(&wl, &data).unwrap();
    // println!("{:?}", spd);
        assert_ulps_eq!(spd[0], 0.);
        assert_ulps_eq!(spd[50], 0.5);
        assert_ulps_eq!(spd[100], 1.0);
        assert_ulps_eq!(spd[250], 0.5);
        assert_ulps_eq!(spd[400], 0.0);

        // top hat with slanted angles
        data = vec![0.0, 1.0, 1.0, 0.0];
        wl = vec![480.0, 490.0, 570.0, 580.0];
        spd = linterp_irr(&wl, &data).unwrap();
    // println!("{:?}", spd);
        assert_ulps_eq!(spd[0], 0.0);
        assert_ulps_eq!(spd[100], 0.0);
        assert_ulps_eq!(spd[110], 1.0);
        assert_ulps_eq!(spd[190], 1.0);
        assert_ulps_eq!(spd[200], 0.0);
        assert_ulps_eq!(spd[300], 0.0);
        assert_ulps_eq!(spd[400], 0.0);

    }

}