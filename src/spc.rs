use std::{collections::BTreeMap, error::Error, ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign}};

use url::Url;
use wasm_bindgen::{prelude::wasm_bindgen, JsValue};

use nalgebra::{DVector, SVector};

use crate::{cri::CRI, data::{A, D50, D65}, obs::Observer, physics::{gaussian_peak_one, led_ohno, planck, stefan_boltzmann}, CmError};


#[wasm_bindgen]
#[derive(PartialEq, Eq, Clone, Copy)]
pub enum Category { 
    /// The spectral distribution of onne or more sources, illuminating a color sample
    Illuminant, 

    /// A Filter spectrum , such as a wratten or glass filter, which changes the properties of an illuminant.
    Filter,     

    /// The spectrum of a color patch, typically consisting of a paint or ink on a substrate, as measured with a spectrophotomteer.
    ColorPatch,

    /// A ray of light from object we are looking at, typically an illuminated by an illuminant.
    Stimulus,   
    
    /// The type of spectrum is unknown.
    Unknown,   
}

// Standard Spectrum domain ranging from 380 to 780 nanometer,
// with 401 values.
pub const NS: usize = 401;

#[derive(Clone, Copy)]
pub enum StdIlluminant {
    D65,
    D50,
    A
}

impl StdIlluminant {
    pub fn spectrum(&self) -> &Spectrum {
        match self {
            StdIlluminant::D65 => &D65,
            StdIlluminant::D50 => &D50,
            StdIlluminant::A => &A,
        }
    }
}

/**
This container holds spectral values within a wavelength domain ranging from 380
to 780 nanometers, with an interval size of 1 nanometer and a total of 401
values. It also includes a category tag and an optional 'total' value for the
aggregate value associated with the spectrum.

The categories are:

- `Illuminant`: a spectral irradiance distribution with values given in watts
    per square meter per nanometer, and a `total` value given in watts per square
    meter.
- `Filter`: a spectral transmission function with unitless values ranging from
    0.0 to 1.0, and the `total` value representing the total power transmission of
    the filter.
- `Substrate`: a spectral transmission function when combined with a `Filter`
    and spectral reflectivity function combined with a `ColorPatch`.
- `ColorPatch`: a spectral reflectivity function with unitless values ranging from
    0.0 to 1.0.
- `Stimulus`: a spectral radiance distribution of a beam of light entering
    through the pupil of our eyes, on its way to be processed and triggering a
    sensation of color in our mind. Spectral data of a stimulus have a unit of watt
    per square meter per nanometer per steradian, and a total.

A `Spectrum` can be constructed from data, but many other construction methods
are available in this library, such as standard illuminants A and D65, Planckian
(Black Body) illuminants, or a `Stimulus` spectrum for a pixel of an sRGB
display.
 */
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub struct Spectrum {
    pub(crate) data: SVector<f64, NS>,
    pub(crate) cat: Category,
    pub(crate) total: Option<f64>, // total irradiance/reflectivity/transimissivity (power based)

}

impl Spectrum {
    pub fn try_new(cat: Category, data: &[f64], total: Option<f64>) -> Result<Self, crate::CmError> {
        if data.len()!=NS {
            Err(crate::CmError::DataSize401Error)
        } else {
            Ok(
                Self {
                    cat,
                    data: SVector::<f64, NS>::from_iterator(data.into_iter().copied()),
                    total
                }
            )
        }
    }
    
    /**
    Get the spectral distribution values as an array.
     */
    pub fn values(&self) -> [f64; NS] {
        self.data.as_slice().try_into().unwrap() // unwrap: we know that data has size NS
    }

    pub fn mul(mut self, rhs: &Self) -> Self {
        self.data = self.data.component_mul(&rhs.data) ;
        self
    }

    pub fn mul_f64(mut self, rhs: f64) -> Self {
        self.data *= rhs;
        self
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

    ```rust
    // Creates a linear gradient filter, with a zero transmission at 380
    // nanometer, and full transmission at 780 nanometer. This is an example
    // using a uniform wavelength domain as input.
    use colorimetry as cmt;
    # use approx::assert_ulps_eq;
    let data = [0.0, 1.0];
    let wl = [380.0, 780.0];
    let mut spd = cmt::Spectrum::linear_interpolate(cmt::Category::Filter, &wl, &data, None).unwrap().values();
    assert_ulps_eq!(spd[0], 0.);
    assert_ulps_eq!(spd[100], 0.25);
    assert_ulps_eq!(spd[200], 0.5);
    assert_ulps_eq!(spd[300], 0.75);
    assert_ulps_eq!(spd[400], 1.0);

    // Creates a top hat filter, with slanted angles, using an irregular
    // wavelength domain.
    let data = vec![0.0, 1.0, 1.0, 0.0];
    let wl = vec![480.0, 490.0, 570.0, 580.0];
    let spd = cmt::Spectrum::linear_interpolate(cmt::Category::Filter, &wl, &data, None).unwrap().values();
    assert_ulps_eq!(spd[0], 0.0);
    assert_ulps_eq!(spd[100], 0.0);
    assert_ulps_eq!(spd[110], 1.0);
    assert_ulps_eq!(spd[190], 1.0);
    assert_ulps_eq!(spd[200], 0.0);
    assert_ulps_eq!(spd[300], 0.0);
    assert_ulps_eq!(spd[400], 0.0);
    ```
    */
    pub fn linear_interpolate(cat: Category, wavelengths: &[f64], data: &[f64], total: Option<f64>) ->Result<Self, crate::CmError> {
        let spdata = match wavelengths.len() {
           2 =>  linterp(wavelengths.try_into().unwrap(), data)?,
           3.. => linterp_irr(wavelengths, data)?,
           _ => return Err(crate::CmError::InterpolateWavelengthError)
        };
        Ok(Self::try_new(cat, &spdata, total)?)
    }

    /**
    Smooth a Spectrum by convolution with a Gaussian function
     */
    pub fn smoothing_filter(mut self, mut width: f64) -> Self {
        if width < 1E-3 { width *= 1E6 }; // to nanometer
        let sigma = width / ((8.0 * 2f64.ln()).sqrt());
        let sd3 = (3.0 * sigma).floor() as i32;
        let kernel =  DVector::<f64>::from_iterator((2*sd3+1) as usize, (-sd3..=sd3).into_iter().map(|i| gaussian_peak_one(i as f64, 0.0, sigma)));
        self.data = self.data.convolve_same(kernel);
        self
    }

    /**
     Standard Daylight Spectrum representing average daylight.

     It's truncated from the official standard, which
     covers 300 to 830 nanometers. It has a correlated color temperature of 6500
     K and should be used in color calculations requiring representative
     daylight. Variations occur based on factors like season, time of day, and
     location. For more details, refer to ISO 10526:1999/CIE
     S005/E-1998.
     */
    pub fn d65_illuminant() -> Self {
        D65.clone()
    }

    /// CIE D50 Illuminant Standard Spectrum with 401 values over a range from
    /// 380 to 780 nanometers, with an interval size of 1 nanometer. Please be
    /// aware that this spectrum is truncated from the official standard, which
    /// is defined over a range from 300 to 830 nanometer.
    ///
    /// For most applications CIE recommends to use the D65 illuminant, to
    /// represent daylight, but this illuminant is often used in the printing
    /// industry.
    pub fn d50_illuminant() -> Self {
        D50.clone()
    }

    /// CIE A Illuminant Standard Spectrum with 401 values over a range from 380 to 780
    /// nanometers, with an interval size of 1 nanometer. This illuminant 
    /// This illuminant is intended to represent typical, domestic,
    /// tungsten-filament lighting. Its relative spectral power distribution is
    /// that of a Planckian radiator at a temperature of approximately 2856 K.
    /// CIE standard illuminant A should be used in all applications of
    /// colorimetry involving the use of incandescent lighting, unless there are
    /// specific reasons for using a different illuminant.
    /// See _— ISO 10526:1999/CIE S005/E-1998, CIE Standard Illuminants for Colorimetry_
    pub fn a_illuminant() -> Self {
        A.clone()
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
    pub fn band_filter(center: f64, width: f64) -> Self {
        let [center_m, width_m] = to_meter([center, width]);
        let data = SVector::<f64,NS>::from_fn(|i,_j|
            {
                let w = (i+380) as f64 * 1E-9;
                let left = center_m - width_m/2.0;
                let right = center_m + width_m/2.0;
                if w < left - f64::EPSILON || w > right + f64::EPSILON { 0.0}
                else {1.0}
            }
        );
        Self { data, cat: Category::Filter, total: None }
    }

    /// A Gaussian Filter, specified by a central wavelength, and a
    /// full-width-half-maximum value, both in units of meter, or nanometer.
    ///
    /// The filter has a peak value of 1.0
    pub fn gaussian_filter(center: f64, width: f64) -> Self {
        let [center_m, width_m] = to_meter([center, width]);
        let data = SVector::<f64,NS>::from_fn(|i,_j|
            gaussian_peak_one((i+380) as f64 * 1E-9, center_m, width_m)
        );
        Self { data, cat: Category::Filter, total: None }
    }

    /// Spectral transmission of a theoretical grey filter, with a constant
    /// transmission over a range from 380 to 780 nanometer, with 1 nanometer
    /// intervals. The transmission value is given as its argument, and should
    /// be in the range from 0.0 to 1.0.  Values outside this range are clamped
    /// to 0.0 for negative values, and 1.0 for values greater than 1.0.
    pub fn gray_filter(gval: f64) -> Self {
        Self{ 
            data: SVector::<f64,NS>::repeat(gval.clamp(0.0, 1.0)),
            cat: Category::Filter,
            total: None,
        }
    }

    /// Theoretical spectrum of a perfect grey color patch, consisting of 401
    /// values equal to the value given in the argument, over a range from 380
    /// to 780 nanometer. Mainly used for color mixing calculations.
    pub fn gray(gval: f64) -> Self {
        Self{ 
            data: SVector::<f64,NS>::repeat(gval.clamp(0.0, 1.0)),
            cat: Category::ColorPatch,
            total: None,
        }
    }

    /// Theoretical spectrum of a perfect white color patch, consisting of 401
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

    /// E, or Equal Energy Illuminant with an irradiance of 1 Watt per square
    /// meter in the spectrum between 380 and 780 nanometer
    pub fn equal_energy_illuminant() -> Self {
        let s = 1./NS as f64;
        Self{ 
            data: SVector::<f64,NS>::repeat(s),
            cat: Category::Illuminant,
            total: Some(1.0)
        }


    }

    /// A pure thermal emission based illuminant according to Planck's law.
    /// 
    /// The generated spectrum is scaled to have a total power, over the full
    /// spectrum (including infrared), of 1 Watt.
    /// ```rust
    /// # use crate::colorimetry::{Spectrum, CIE1931};
    /// # use approx::assert_ulps_eq;
    /// 
    /// let p3000 = Spectrum::planckian_illuminant(3000.0);
    /// let [x, y] = CIE1931.xyz(&p3000).chromaticity();
    /// let l = CIE1931.xyz(&p3000).luminous_value();
    /// assert_ulps_eq!(l, 20.668_927, epsilon = 1E-6);
    /// assert_ulps_eq!(x, 0.436_935, epsilon = 1E-6);
    /// assert_ulps_eq!(y, 0.404_083, epsilon = 1E-6);
    /// 
    /// ```
    /// ```javascript
    /// let x = 350.0;
    /// ```
    pub fn planckian_illuminant(cct: f64) -> Self {

        let s = 1E-9/stefan_boltzmann(cct); // 1W/m2 total irradiance
        let data = SVector::<f64,NS>::from_fn(|i,_j|s * planck((i+380) as f64*1e-9, cct));
        Self {
            data,
            cat: Category::Illuminant,
            total: Some(1.0),
        }
    }


    /// A spectral power distribution for a Light Emitting Diode.
    ///
    /// The spectrum is definded by a center wavelength, in units of meter or
    /// nanometer, and a full-width-half-maximum value, also in units of meter
    /// or nanometer. The generated spectrum is based on the model as published
    /// by Yoshi Ohno, from NIST, in his article, _Spectral Design
    /// considerations for white LED Color Rendering_, **Optical Engineering 44(11)**, 
    /// November 2005.
    pub fn led_illuminant(center: f64, width: f64) -> Self {
        let [center_m, width_m] = to_meter([center, width]);
        let data = SVector::<f64,NS>::from_fn(|_i,j|
            led_ohno((j+380) as f64 * 1E-9, center_m, width_m) * 1E-9
        );
        Self { data, cat: Category::Illuminant, total: Some(1.0) }
    }

    // Sets irradiance, tyically expressed in units of Watt per square meter.
    // Also overwrite spectrum type to irradiance.
    pub fn set_irradiance(&mut self, irradiance: f64) -> &mut Self {
        let s = if let Some(t) = self.total {
            irradiance/t
        } else {
            irradiance/self.data.sum()
        };
        self.data.iter_mut().for_each(|v|*v = *v *s);
        self.cat = Category::Illuminant;
        self.total = Some(irradiance);
        self
    }

    // Calculate a spectrum's irradiance if it is an illuminant.
    // Produces a "Not A Number" value, if not an illuminant.
    pub fn irradiance(&self) -> f64 {
        if self.cat == Category::Illuminant && self.total.is_some() {
            self.total.unwrap()
        } else {
            f64::NAN
        }
    }

    pub fn set_illuminance(mut self, obs: &Observer, illuminance: f64) -> Self {
        let l = illuminance / (obs.data.row(1) *  self.data * obs.lumconst).x;
        self.data.iter_mut().for_each(|v| *v = *v * l);
        self.cat = Category::Illuminant;
        self
    }

    pub fn illuminance(&self, obs: &Observer) -> f64 {
        if self.cat == Category::Illuminant {
            (obs.data.row(1) * self.data *  obs.lumconst).x
        } else {
            f64::NAN
        }
    }

    /// Downloads a spectrum
    pub async fn fetch(loc: &str) -> Result<Self, Box<dyn Error>> {
        let _url = Url::parse(loc)?;
        todo!()


    }

    /// Calculates the Color Rendering Index values for illuminant spectrum.
    /// 
    /// To use this function, first use `CRI::init().await`, which downloads the
    /// Test Color Samples required for the calculation.  These are downloaded
    /// seperately to limit the size of the main web assembly library.
    pub fn cri(&self) -> Result<CRI, Box<dyn Error>> {
        todo!()
    }
}


#[cfg(test)]
mod spectrum_tests {
    use crate::{Spectrum, D65, CIE1931};
    use approx::assert_ulps_eq;

    #[test]
    fn test_led(){
        use approx::assert_ulps_eq;
        let ls = Spectrum::led_illuminant(550.0, 25.0);
        assert_ulps_eq!(ls.irradiance(), 1.0, epsilon = 1E-9);
    }

    #[test]
    fn test_chromaticity(){
        let xyz0 = CIE1931.xyz(&D65);
        let [x0, y0] = xyz0.chromaticity();

        let illuminance = D65.illuminance(&CIE1931);
        let d65 = D65.set_illuminance(&CIE1931, 100.0);
        let xyz = CIE1931.xyz(&d65);
        let [x, y] = xyz.chromaticity();
    
        assert_ulps_eq!(x0, x);
        assert_ulps_eq!(y0, y);
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
    pub fn new_js(cat: Category, data: &[f64], total: Option<f64>) -> Result<Spectrum, wasm_bindgen::JsError> {
        Ok(Spectrum::try_new(cat, data, total)?)
    }

    /// Returns the spectral data values, as a Float64Array containing 401 data
    /// points, over a wavelength domain from 380 t0 780 nanometer, with a
    /// stepsize of 1 nanometer.
    #[wasm_bindgen(js_name=Values)]
    pub fn values_js(&self) -> Box<[f64]> {
        self.values().into()
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
    let mut spd = cmt::Spectrum::linear_interpolate(cmt::Category::Filter, &wl, &data, None).unwrap().values();
    assert_ulps_eq!(spd[0], 0.);
    assert_ulps_eq!(spd[100], 0.25);
    assert_ulps_eq!(spd[200], 0.5);
    assert_ulps_eq!(spd[300], 0.75);
    assert_ulps_eq!(spd[400], 1.0);

    // Creates a top hat filter, with slanted angles, using an irregular
    // wavelength domain.
    let data = vec![0.0, 1.0, 1.0, 0.0];
    let wl = vec![480.0, 490.0, 570.0, 580.0];
    let spd = cmt::Spectrum::linear_interpolate(cmt::Category::Filter, &wl, &data, None).unwrap().values();
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
    pub fn linear_interpolate_js(cat: Category, wavelengths: &[f64], data: &[f64], total_js: &JsValue) -> Result<Spectrum, crate::CmError> {
        let total = if let Some(t) = total_js.as_f64() {
            Some(t)
        } else {
            None
        };
        Self::linear_interpolate(cat, wavelengths, data, total)

    }

    /// Calculates the Color Rendering Index values for illuminant spectrum.
    /// 
    /// To use this function, first use `await CRI.init()`, which downloads the
    /// Test Color Samples required for the calculation.  These are downloaded
    /// seperately to limit the size of the main web assembly library.
    #[wasm_bindgen(js_name=cri)]
    pub fn cri_js(&self) -> Result<CRI, crate::CmError> {
        todo!()
    }

}

fn mixed_category(s1: &Spectrum, s2: &Spectrum) -> Category {
    if 
        (s1.cat == Category::Illuminant && s2.cat == Category::ColorPatch) ||
        (s1.cat == Category::ColorPatch && s2.cat == Category::Illuminant) ||
        (s1.cat == Category::Illuminant && s2.cat == Category::Filter) ||
        (s1.cat == Category::Filter && s2.cat == Category::Illuminant) {
        return Category::Stimulus
    } else if s1.cat == Category::Filter && s2.cat == Category::Filter {
            Category::Filter
    } else if s1.cat == Category::ColorPatch && s2.cat == Category::ColorPatch {
            Category::ColorPatch
    } else {
            Category::Unknown
    }
}

// Multiplication of two spectra, typically for a combinations of an illuminant and a filter or ColorPatch,
// or when combining multiple ColorPatchs or filters. Subtractive Mixing.
impl Mul for Spectrum {
    type Output = Self;

    // multiply two cie spectra
    fn mul(self, rhs: Self) -> Self::Output {
        Self{
            cat: mixed_category(&self, &rhs), 
            data: self.data.component_mul(&(rhs.data)),
            total: None,
        }
    }
}

impl Mul<f64> for Spectrum {
    /// Multiply a spectrum with a scalar f64 value.
    /// ```
    ///     use crate::colorimetry::Spectrum;
    ///     use approx::assert_ulps_eq;
    ///
    ///     let mut led = Spectrum::led_illuminant(550.0, 25.0);
    ///     let mut irradiance = led.irradiance();
    ///     assert_ulps_eq!(led.irradiance(), 1.0, epsilon = 1E-10);
    ///
    ///     led = led * 10.0;
    ///     assert_ulps_eq!(led.irradiance(), 10.0, epsilon = 1E-10);
    /// ```
    type Output = Spectrum;

    // spectrum * scalar
    fn mul(self, rhs: f64) -> Self::Output {
        let total = if let Some(t) = self.total {
            Some(t * rhs)
        } else {
            None
        };
        Self{
            cat: self.cat,
            data: self.data * rhs,
            total
        }
    }
}

impl Mul<Spectrum> for f64 {
    /// Multiply a spectrum with a scalar f64 value.
    /// ```
    ///     use crate::colorimetry::Spectrum;
    ///     use approx::assert_ulps_eq;
    ///
    ///     let mut led = Spectrum::led_illuminant(550.0, 25.0);
    ///     let mut irradiance = led.irradiance();
    ///     assert_ulps_eq!(led.irradiance(), 1.0, epsilon = 1E-10);
    ///
    ///     led = 10.0 * led;
    ///     assert_ulps_eq!(led.irradiance(), 10.0, epsilon = 1E-10);
    /// ```
    type Output = Spectrum;

    // scalar * spectrum
    fn mul(self, rhs: Spectrum) -> Self::Output {
        let total = if let Some(t) = rhs.total {
            Some(t * self)
        } else {
            None
        };
        Self::Output {
            cat: rhs.cat,
            data: self * rhs.data,
            total
        }
    }
}

#[test]
fn mul_spectra_test(){
    use approx::assert_ulps_eq;
    let g = Spectrum::gray(0.5);

    let w = 2.0 * g;
    for i in 380..780 {
        assert_ulps_eq!(w[i], 1.0);
    }

    let v = g * 2.0;
    for i in 380..780 {
        assert_ulps_eq!(v[i], 1.0);
    }

}


// Addition of spectra, typically used for illuminant (multiple sources).
// Additive mixing
impl Add for Spectrum {
    type Output = Self;

    // add two cie spectra
    fn add(self, rhs: Self) -> Self::Output {
        let category = if self.cat == Category::Illuminant && rhs.cat == Category::Illuminant {
            Category::Illuminant
        } else {
            Category::Unknown
        };
        Self{
            cat: category,
            data: self.data + rhs.data,
            total: None
        }
    }
}

// Addition of spectra, typically used for illuminant (multiple sources).
// Additive mixing
impl AddAssign for Spectrum {
    fn add_assign(&mut self, rhs: Self) {
        self.data += rhs.data
    }
}

#[test]
fn add_spectra(){
    use approx::assert_ulps_eq;
    let mut g1 = Spectrum::gray(0.5);
    let g2 = Spectrum::gray(0.5);
    let g = g1 + g2;
    for i in 380..780 {
        assert_ulps_eq!(g[i], 1.0);
    }

    g1 += g2;
    for i in 380..780 {
        assert_ulps_eq!(g1[i], 1.0);
    }

    let v = 2.0 * Spectrum::gaussian_filter(550.0, 50.0) + -2.0 * Spectrum::gaussian_filter(550.0, 50.0);
    for i in 380..780 {
        assert_ulps_eq!(v[i], 0.0);
    }

}

impl MulAssign for Spectrum {
    /// Element wise multiply (filter) a spectrum with another spectrum.
    fn mul_assign(&mut self, rhs: Self) {
        self.data.iter_mut().zip(rhs.data.iter()).for_each(|(v,w)| *v *= *w);

    }
}

impl MulAssign<f64> for Spectrum {
    /// Scale a spectrum with a scaler value.
    /// Depending on the type of spectrum this has different meanings.
    /// - for an illuminant, this scales the irradiance,
    /// - for a ColorPatch, this scales the total reflectivity.
    /// - for a filter, it changes its transmission.
    /// ```
    ///     use crate::colorimetry::Spectrum;
    ///     use approx::assert_ulps_eq;
    ///
    ///     let mut led = Spectrum::led_illuminant(550.0, 25.0);
    ///     assert_ulps_eq!(led.irradiance(), 1.0, epsilon = 1E-10);
    ///
    ///     led *= 10.0;
    ///     assert_ulps_eq!(led.irradiance(), 10.0, epsilon = 1E-10);
    /// ```
    fn mul_assign(&mut self, rhs: f64) {
        self.total = if let Some(t) = self.total {
            Some(t * rhs)
        } else {
            None
        };
        self.data.iter_mut().for_each(|v| *v *= rhs);

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
            &self.data[(i-380,0)]

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
            380..=780 => &mut self.data[(i-380,0)],
            ..380 => &mut self.data[(0,0)],
            781.. => &mut self.data[(400,0)],
        }
    }
}
#[test]
fn index_test(){
    use approx::assert_ulps_eq;
    let mut s = Spectrum::white();

    // Set a spectral value
    s[500] = 0.5;
    assert_ulps_eq!(s[500], 0.5);

    // A read from an index out of the 380..=780 range with given f64::NAN.
    s[300] = 0.5;
    assert!(s[300].is_nan());
}

/// Convenience function for specifying wavelengths in nanometers or meters.
/// Wavelength values larger than !E-3 are assumed to have the unit nanometer
/// and are converted to a unit of meters.
fn to_meter<const N: usize>(mut v:[f64; N]) -> [f64;N] {
    if v.iter().any(|v|v>&1E-3) {
        v.iter_mut().for_each(|v| *v *= 1E-9)
    };
    v
}



#[test]
fn test_to_meter() {
    use approx::assert_ulps_eq;

    let mut v1 = [380.0];
    v1 = to_meter(v1);
    assert_ulps_eq!(v1[0], 380E-9);

    let mut v2 = [380E-9, 780E-9];
    v2 = to_meter(v2);
    assert_ulps_eq!(v2[0], 380E-9);
    assert_ulps_eq!(v2[1], 780E-9);

}



#[cfg(test)]
mod tests {

    use crate::spc::Spectrum;

    use crate::data::CIE1931;
    use approx::assert_ulps_eq;

    #[test]
    fn ee() {
        let [x, y ] = CIE1931.xyz(&&Spectrum::equal_energy_illuminant().set_illuminance(&CIE1931, 100.0)).chromaticity();
        assert_ulps_eq!(x, 0.333_3, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.333_3, epsilon = 5E-5);
    }

    #[test]
    fn d65() {
        let [x, y ] = CIE1931.xyz(&Spectrum::d65_illuminant().set_illuminance(&CIE1931, 100.0)).chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(x, 0.312_72, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.329_03, epsilon = 5E-5);
    }

    #[test]
    fn d50() {
        let [x, y ] = CIE1931.xyz(&Spectrum::d50_illuminant().set_illuminance(&CIE1931, 100.0)).chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(x, 0.345_67, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.358_51, epsilon = 5E-5);
    }

    #[test]
    fn a() {
        let [x, y ] = CIE1931.xyz(&Spectrum::a_illuminant().set_illuminance(&CIE1931, 100.0)).chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(x, 0.447_58, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.407_45, epsilon = 5E-5);
    }
}

/// Linear interpolatino over a dataset over an equidistant wavelength domain
fn linterp(mut wl: [f64;2], data: &[f64]) -> Result<[f64;NS], CmError> {
    wl.sort_by(|a, b| a.partial_cmp(b).unwrap()); 
    let [wl, wh] = to_meter(wl);
    let dlm1 = data.len()-1; // data length min one
    
    let mut spd = [0f64; NS];
    spd.iter_mut().enumerate().for_each(|(i,v)|{
        let l = (i + 380) as f64 * 1E-9;
        let t = ((l-wl)/(wh - wl)).clamp(0.0, 1.0);
        let j = (t * dlm1 as f64).trunc() as usize;
        let f = t.fract();
        if j >= dlm1 {
            *v = data[dlm1];
        } else {
            *v = data[j] * (1.0 - f) + data[j+1] * f;
        }
    });
    Ok(spd)
}

#[test]
fn test_linterp(){
    use approx::assert_ulps_eq;

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

/**
Spectrum constructed by linear interpolatino over a dataset with an irregular
wavelength domain.

This algorithm uses a BTreeMap coolection, with wavelengths in picometers as key,
to find a data interval containing the target wavelengths.
 */
fn linterp_irr(wl: &[f64], data: &[f64]) -> Result<[f64;NS], CmError> {
    if wl.len()!=data.len() {
        Err(CmError::InterpolateWavelengthError)
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
