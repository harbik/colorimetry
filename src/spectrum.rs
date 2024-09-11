/*!
# Spectral Composition

The field of Colorimetry uses mathematical models to describe the sensations in our mind which we
call color.  These models are based the spectral composition of stimuli, essentially rays of light
hitting the photosensitive cells in the back of our eyes, and the spectral sensitiviy of these 
cells. The spectral composition of light, and the objects involved in its processing such as filters
and painted patches, is represented by the [Spectrum]-object in this library.
The spectral sensitivity of human vision is described by an [Observer](crate::Observer).
*/
use std::{collections::BTreeMap, default, error::Error, iter::Sum, ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign}};

use num_traits::ToPrimitive;
use url::Url;
use wasm_bindgen::{prelude::wasm_bindgen, JsValue};

use nalgebra::{DVector, SVector};

use crate::{data::{D50, D65}, observer::ObserverData, physics::{gaussian_peak_one, led_ohno, planck, stefan_boltzmann}, wavelength, CmtError, StdIlluminant, C, CIE1931, RGB};



/**
A `Category` is used as a tag in a `Spectrum`, for use for operations specific for a particular
spectrum type, and to avoid incorrect colorimetric calculations.

- `Illuminant`: a spectral irradiance distribution with values given in watts
    per square meter per nanometer, and a `total` value given in watts per square
    meter.

- `Colorant`: a spectral reflectivity (color patch) or transmissivity with unitless values ranging from
    0.0 to 1.0, which changes an illuminant's spectral composition to produce a stimulus.

- `Stimulus`: a spectral radiance distribution of a beam of light entering
    through the pupil of our eyes, on its way to be processed and triggering a
    sensation of color in our mind. Spectral data of a stimulus have a unit of watt
    per square meter per nanometer per steradian, and a total.

 */
#[wasm_bindgen]
#[derive(PartialEq, Eq, Clone, Copy, Default, Debug)]
pub enum Category { 
    /// The spectral distribution of onne or more sources, illuminating a color sample
    Illuminant, 

    /// The spectrum of a color patch, typically consisting of a paint or ink on a substrate, as measured with a spectrophotomteer.
    Colorant,

    /// A ray of light from object we are looking at, typically an illuminated by an illuminant.
    Stimulus,   
    
    /// The type of spectrum is unknown.
    #[default]
    Unknown,   
}

// Standard Spectrum domain ranging from 380 to 780 nanometer,
// with 401 values.
pub const NS: usize = 401;

/// Multipilication definition as used in the Spectum multiplication operator.
///
/// Multiplication of same categories results in same.
/// Illuminants with filters or patches result in stimuli.
impl Mul for Category {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self == rhs {
            return self
        } else {
            match (self, rhs) {
                (Category::Illuminant, Category::Colorant) => Category::Stimulus,
                (Category::Colorant, Category::Illuminant) => Category::Stimulus,
                (Category::Stimulus, _ ) => Category::Stimulus,
                (_, Category::Stimulus) => Category::Stimulus,
                _ => Category::Unknown,
            }
        }
    }
}

/// Addition definition as used in the Spectum multiplication operator.
///
/// Addition of same categories stay the same.  A combination of a filter and a patch results in a
/// "filtered patch', so a patch... If you add filters and/or patches, Anything else is undefined.
impl Add for Category {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self == rhs {
            self
        } else { 
            // can not add different categories
            Category::Unknown
        }
    }
}

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
pub struct Spectrum {
    pub(crate) data: SVector<f64, NS>,
    pub(crate) cat: Category,

}

impl Spectrum {
    pub fn try_new(cat: Category, data: &[f64]) -> Result<Self, crate::CmtError> {
        if data.len()!=NS {
            Err(crate::CmtError::DataSize401Error)
        } else {
            if cat == Category::Colorant {
                Self::try_new_colorant(data)
            } else {
                Ok(
                    Self {
                        cat,
                        data: SVector::<f64, NS>::from_iterator(data.into_iter().copied()),
                    }
                )
            }
        }
    }

    /// Create a Colorant Spectrum, with data check.
    /// 
    /// In this library, the Colorant category represents color filters and color patches, with have spectral values between 0.0 and 1.0.
    /// This is the preferred way to crete a colorant, as it does a range check.
    pub fn try_new_colorant(data: &[f64]) -> Result<Self, crate::CmtError> {
        if data.iter().any(|&v|v<0.0 || v>1.0){
            Err(CmtError::OutOfRange { name: "Colorant Spectral Value".into(), low: 0.0, high: 1.0 })
        } else {
                Ok(
                    Self {
                        cat: Category::Colorant,
                        data: SVector::<f64, NS>::from_iterator(data.into_iter().copied()),
                    }
                )
        }
    }
    
    pub fn try_set_category(mut self, cat: Category) -> Result<Self, CmtError> {
        if self.data.iter().any(|&v|v<0.0 || v>1.0){
            Err(CmtError::OutOfRange { name: "Colorant Spectral Value".into(), low: 0.0, high: 1.0 })
        } else {
            self.cat = cat;
            Ok(self)
        }
    }

    pub fn set_category_unchecked(mut self, cat: Category) -> Self {
        self.cat = cat;
        self
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
    let mut spd = cmt::Spectrum::linear_interpolate(cmt::Category::Colorant, &wl, &data).unwrap().values();
    assert_ulps_eq!(spd[0], 0.);
    assert_ulps_eq!(spd[100], 0.25);
    assert_ulps_eq!(spd[200], 0.5);
    assert_ulps_eq!(spd[300], 0.75);
    assert_ulps_eq!(spd[400], 1.0);

    // Creates a top hat filter, with slanted angles, using an irregular
    // wavelength domain.
    let data = vec![0.0, 1.0, 1.0, 0.0];
    let wl = vec![480.0, 490.0, 570.0, 580.0];
    let spd = cmt::Spectrum::linear_interpolate(cmt::Category::Colorant, &wl, &data).unwrap().values();
    assert_ulps_eq!(spd[0], 0.0);
    assert_ulps_eq!(spd[100], 0.0);
    assert_ulps_eq!(spd[110], 1.0);
    assert_ulps_eq!(spd[190], 1.0);
    assert_ulps_eq!(spd[200], 0.0);
    assert_ulps_eq!(spd[300], 0.0);
    assert_ulps_eq!(spd[400], 0.0);
    ```
    */
    pub fn linear_interpolate(cat: Category, wavelengths: &[f64], data: &[f64]) ->Result<Self, crate::CmtError> {
        let spdata = match wavelengths.len() {
           2 =>  linterp(wavelengths.try_into().unwrap(), data)?,
           3.. => linterp_irr(wavelengths, data)?,
           _ => return Err(crate::CmtError::InterpolateWavelengthError)
        };
        Ok(Self::try_new(cat, &spdata)?)
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
        Self { data, cat: Category::Colorant}
    }

    /// A Gaussian Filter, specified by a central wavelength, and a
    /// full-width-half-maximum value, both in units of meter, or nanometer.
    ///
    /// The filter has a peak value of 1.0
    pub fn gaussian_filter(center: f64, width: f64) -> Self {
        let [center_m, width_m] = wavelengths([center, width]);
        let data = SVector::<f64,NS>::from_fn(|i,_j|
            gaussian_peak_one((i+380) as f64 * 1E-9, center_m, width_m)
        );
        Self { data, cat: Category::Colorant }
    }

    /// Theoretical spectrum of a perfect grey colorant, consisting of 401
    /// values equal to the value given in the argument, over a range from 380
    /// to 780 nanometer. Mainly used for color mixing calculations.
    pub fn gray(gval: f64) -> Self {
        Self{ 
            data: SVector::<f64,NS>::repeat(gval.clamp(0.0, 1.0)),
            cat: Category::Colorant,
        }
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


    /// A spectral composition of a display pixel, set to three sRGB color values.  The spectrum is
    /// a linear combination of the spectral primaries, which are Gaudssian filtered components in
    /// this library.
    pub fn srgb(r_u8: u8, g_u8: u8, b_u8: u8) -> Self {
        let rgb = RGB::from_u8(r_u8, g_u8, b_u8, Some(crate::Observer::Std1931), Some(crate::RgbSpace::SRGB));
        rgb.into()
    }

    /// A spectral composition of a display pixel, set to three sRGB color values.  The spectrum is
    /// a linear combination of the spectral primaries, which are Gaudssian filtered components in
    /// this library.
    pub fn rgb(rgb: RGB) -> Self {
        rgb.into()
    }

    /// E, or Equal Energy Illuminant with an irradiance of 1 Watt per square
    /// meter in the spectrum between 380 and 780 nanometer
    pub fn equal_energy_illuminant() -> Self {
        let s = 1./NS as f64;
        Self{ 
            data: SVector::<f64,NS>::repeat(s),
            cat: Category::Illuminant,
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
        let [center_m, width_m] = wavelengths([center, width]);
        let data = SVector::<f64,NS>::from_fn(|i,_j|
            // j = 0, first column
            led_ohno(wavelength(i+380), center_m, width_m) * 1E-9
        );
        Self { data, cat: Category::Illuminant }
    }

    // Sets irradiance, tyically expressed in units of Watt per square meter.
    // Also overwrite spectrum type to Illuminant
    pub fn set_irradiance(mut self, irradiance: f64) -> Self {
        let s = irradiance/self.data.sum();
        self.data.iter_mut().for_each(|v|*v = *v *s);
        self.cat = Category::Illuminant;
        self
    }

    // Calculate a spectrum's irradiance if it is an illuminant.
    // Produces a "Not A Number" value, if not an illuminant.
    pub fn irradiance(&self) -> f64 {
        if self.cat == Category::Illuminant {
            self.data.sum()
        } else {
            f64::NAN
        }
    }

    pub fn set_illuminance(mut self, obs: &ObserverData, illuminance: f64) -> Self {
        let l = illuminance / (obs.data.row(1) *  self.data * obs.lumconst).x;
        self.data.iter_mut().for_each(|v| *v = *v * l);
        self.cat = Category::Illuminant;
        self
    }

    pub fn illuminance(&self, obs: &ObserverData) -> f64 {
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
    #[cfg(feature="cri")]
    pub fn cri(&self) -> Result<crate::CRI, CmtError> {
        self.try_into()
    }

}

impl Default for Spectrum {
    fn default() -> Self {
        Self { 
            data: SVector::<f64, NS>::zeros(),
            cat: Default::default(), 
         }
    }
}

impl Sum for Spectrum {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut s = Self::default() ;
        iter.for_each(|si|s += si);
        s
    }
}

/// Spectral representation the color of a display pixel, described by a [RGB]
/// instance.
///
/// It uses a linear combination of the spectral primaries as defined for a particular
/// [RgbSpace](crate::RgbSpace).
/// Most of the color spaces in this library use Daylight filtered Gaussian primaries,
/// but you can also use your own color space based on primaries measured by a spectrometer.
/// Spectral representations of pixels allow color matching for arbitrary observers,
/// not only the CIE 1931 standard observer.
impl From<RGB> for Spectrum {
    fn from(rgb: RGB) -> Self {
        let prim = &rgb.space.data().0.primaries;
        let yrgb = rgb.observer.data().rgb2xyz(&rgb.space).row(1);
        let mut s: Spectrum = rgb.data.iter().zip(yrgb.iter()).zip(prim.iter()).map(|((v,w),s)|*v * *w * *s).sum();
        s.cat = Category::Stimulus;
        s
    }
}

#[cfg(test)]
mod spectrum_tests {
    use crate::{Spectrum, D65, CIE1931, RGB};
    use approx::assert_ulps_eq;

    #[test]
    fn test_spectrum_from_rgb(){
        let white = RGB::new(1.0, 1.0, 1.0, None, None).into();
        approx::assert_ulps_eq!(CIE1931.xyz(&white, None), CIE1931.xyz_d65().set_illuminance(100.0), epsilon = 1E-6);
        let red = Spectrum::srgb(255, 0, 0);
        assert_ulps_eq!(CIE1931.xyz(&red, None).chromaticity().as_ref(), &[0.64, 0.33].as_ref(), epsilon = 1E-5);
    }

    #[test]
    fn test_led(){
        use approx::assert_ulps_eq;
        let ls = Spectrum::led_illuminant(550.0, 25.0);
        assert_ulps_eq!(ls.irradiance(), 1.0, epsilon = 1E-9);
    }

    #[test]
    fn test_chromaticity(){
        let xyz0 = CIE1931.xyz(&D65, None);
        let [x0, y0] = xyz0.chromaticity();

        let illuminance = D65.illuminance(&CIE1931);
        let d65 = D65.set_illuminance(&CIE1931, 100.0);
        let xyz = CIE1931.xyz(&d65, None);
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
    pub fn new_js(cat: Category, data: &[f64]) -> Result<Spectrum, wasm_bindgen::JsError> {
        Ok(Spectrum::try_new(cat, data)?)
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
    pub fn linear_interpolate_js(cat: Category, wavelengths: &[f64], data: &[f64], total_js: &JsValue) -> Result<Spectrum, crate::CmtError> {
        Self::linear_interpolate(cat, wavelengths, data)

    }

    /// Calculates the Color Rendering Index values for illuminant spectrum.
    /// 
    /// To use this function, first use `await CRI.init()`, which downloads the
    /// Test Color Samples required for the calculation.  These are downloaded
    /// seperately to limit the size of the main web assembly library.
    #[cfg(feature="cri")]
    #[wasm_bindgen(js_name=cri)]
    pub fn cri_js(&self) -> Result<crate::cri::CRI, crate::CmtError> {
        todo!()
    }

    /// Get the StdIlluminant spectrum. Typically you don't need to use the Spectrum itself, as many
    /// methods just accept the StdIlluminant directly.
    #[wasm_bindgen(js_name=illuminant)]
    pub fn llluminant_js(stdill: StdIlluminant) -> Self {
        // need this as wasm_bindgen does not support `impl` on Enum types (yet?).
        // in Rust use StdIlluminant.spectrum() directly, which also gives a reference instead of a copy.
        *stdill.spectrum()
    }
}


// Multiplication of two spectra using the `*`-operator, typically for a combinations of an illuminant and a colorant
// or when combining multiple ColorPatchs or filters. Subtractive Mixing.
impl Mul for Spectrum {
    type Output = Self;

    // multiply two cie spectra
    fn mul(self, rhs: Self) -> Self::Output {
        Self{
          //  cat: mixed_category(&self, &rhs), 
            cat: self.cat * rhs.cat,
            data: self.data.component_mul(&(rhs.data)),
        }
    }
}

// Multiplication of two references to spectra using the `*`-operator, typically for a combinations of an illuminant and a colorant,
// or when combining multiple colorants, and producing a new spectrum. Subtractive Mixing.
impl Mul<&Spectrum> for &Spectrum {
    type Output = Spectrum;

    // multiply two cie spectra
    fn mul(self, rhs: &Spectrum) -> Self::Output {
        Self::Output{
          //  cat: mixed_category(&self, &rhs), 
            cat: self.cat * rhs.cat,
            data: self.data.component_mul(&(rhs.data)),
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
        Self{
            cat: self.cat,
            data: self.data * rhs,
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
        Self::Output {
            cat: rhs.cat,
            data: self * rhs.data,
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
        Self{
            cat: self.cat + rhs.cat,
            data: self.data + rhs.data,
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
    /// ```
    /// use colorimetry::{Spectrum, CIE1931, XYZ};
    /// let mut spc = Spectrum::d65_illuminant().set_illuminance(&CIE1931, 100.0);
    /// spc *= Spectrum::white(); // no change in color point, multiply with all ones
    /// let xyz = CIE1931.xyz(&spc); // calculate tristimulus values
    /// approx::assert_ulps_eq!(xyz, CIE1931.xyz_d65());
    /// 
    /// ```
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
/// Can use integer and float values.
/// Wavelength values larger than 1E-3 are assumed to have the unit nanometer
/// and are converted to a unit of meters.
fn wavelengths<T: ToPrimitive, const N: usize>(v:[T; N]) -> [f64;N] {
    v.map(|x|wavelength(x))
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



#[cfg(test)]
mod tests {

    use crate::spectrum::Spectrum;

    use crate::data::CIE1931;
    use crate::StdIlluminant;
    use approx::assert_ulps_eq;

    #[test]
    fn ee() {
        let [x, y ] = CIE1931.xyz(
            &Spectrum::equal_energy_illuminant().set_illuminance(&CIE1931, 100.0), None).chromaticity();
        assert_ulps_eq!(x, 0.333_3, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.333_3, epsilon = 5E-5);
    }

    #[test]
    fn d65() {
        let [x, y ] = CIE1931.xyz(
            &Spectrum::d65_illuminant().set_illuminance(&CIE1931, 100.0), None).chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(x, 0.312_72, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.329_03, epsilon = 5E-5);
    }

    #[test]
    fn d50() {
        let [x, y ] = CIE1931.xyz(&Spectrum::d50_illuminant().set_illuminance(&CIE1931, 100.0), None).chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(x, 0.345_67, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.358_51, epsilon = 5E-5);
    }

    #[cfg_attr(test, cfg(feature="cie-illuminants"))]
    fn a() {
        let [x, y ] = CIE1931.xyz(&StdIlluminant::A.spectrum().set_illuminance(&CIE1931, 100.0), None).chromaticity();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(x, 0.447_58, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.407_45, epsilon = 5E-5);
    }
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
