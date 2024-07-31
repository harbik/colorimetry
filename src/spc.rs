use std::{error::Error, ops::{Add, Mul, MulAssign}};

use url::Url;
use wasm_bindgen::prelude::wasm_bindgen;

use nalgebra::SMatrix;

use crate::{obs::Observer, data::{A, D50, D65}, physics::{led_ohno, planck, stefan_boltzmann}};


#[derive(PartialEq, Eq, Clone, Copy)]
pub enum Category { 
    Illuminant, // one or more sources, illuminating a color sample
    Filter,     // Light Filter, such as a wratten or glass filter
    Substrate,  // e.g. Paper, Glass
    Colorant,   // e.g. Pigment, Paint, Ink, viewed in reflection
    Stimulus,   // A ray of light from object we are looking at, ultimately
                // creating a sensation of color in our mind
    Unknown,   
}

// Standard Spectrum domain ranging from 380 to 780 nanometer,
// with 401 values.
pub const NS: usize = 401;

//pub type StdRowVec = RowVector::<f64, Const<NS>, ArrayStorage<f64, 1, NS>>;
pub type SpDataVec = SMatrix<f64, 1, NS>;

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub struct Spectrum {
    pub(crate) data: SpDataVec,
    pub(crate) cat: Category,
    pub(crate) total: Option<f64>, // total irradiance/reflectivity/transimissivity (power based)

}


#[wasm_bindgen]
impl Spectrum {
    pub fn d65() -> Self {
        D65.clone()
    }

    pub fn d50() -> Self {
        D50.clone()
    }

    pub fn a() -> Self {
        A.clone()
    }

    pub fn grey(gval: f64) -> Self {
        Self{ 
            data: SpDataVec::repeat(gval.clamp(0.0, 1.0)),
            cat: Category::Colorant,
            total: None,
        }
    }

    pub fn grey_filter(gval: f64) -> Self {
        Self{ 
            data: SpDataVec::repeat(gval.clamp(0.0, 1.0)),
            cat: Category::Filter,
            total: None,
        }
    }

    pub fn white() -> Self {
        Self::grey(1.0)
    }

    pub fn black() -> Self {
        Self::grey(0.0)
    }

    // E, or Equal Energy Illuminant with an irradiance of 1 Watt per square
    // meter in the spectrum between 380 and 780 nanometer
    pub fn equal_energy() -> Self {
        let s = 1./NS as f64;
        Self{ 
            data: SpDataVec::repeat(s),
            cat: Category::Illuminant,
            total: Some(1.0)
        }


    }

    pub fn led(center: f64, width: f64) -> Self {
        let [center_m, width_m] = to_meter([center, width]);
        let data = SpDataVec::from_fn(|_i,j|
            led_ohno((j+380) as f64 * 1E-9, center_m, width_m) * 1E-9
        );
        Self { data, cat: Category::Illuminant, total: Some(1.0) }
    }

}

    #[test]
    fn test_led(){
        use approx::assert_ulps_eq;
        let ls = Spectrum::led(550.0, 25.0);
        assert_ulps_eq!(ls.irradiance(), 1.0, epsilon = 1E-9);
    }


impl Spectrum {

    /// A pure thermal emission based illuminant according to Planck's law.
    /// 
    /// The generated spectrum is scaled to have a total power, over the full
    /// spectrum (including infrared), of 1 Watt.
    /// ```
    /// use crate::colorimetry::{Spectrum, CIE1931};
    /// use approx::assert_ulps_eq;
    /// 
    /// let p3000 = Spectrum::planckian(3000.0);
    /// let [l, x, y] = CIE1931.xyz(&p3000).lxy();
    /// assert_ulps_eq!(l, 20.668_927, epsilon = 1E-6);
    /// assert_ulps_eq!(x, 0.436_935, epsilon = 1E-6);
    /// assert_ulps_eq!(y, 0.404_083, epsilon = 1E-6);
    /// 
    /// ```
    pub fn planckian(cct: f64) -> Self {
        let s = 1E-9/stefan_boltzmann(cct); // 1W/m2 total irradiance
        let data = SpDataVec::from_fn(|_i,j|s * planck((j+380) as f64*1e-9, cct));
        Self {
            data,
            cat: Category::Illuminant,
            total: Some(1.0),
        }
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

    pub fn set_illuminance(&mut self, obs: &Observer, illuminance: f64) -> &mut Self {
        let l = illuminance / (self.data * obs.data.column(1) * obs.lumconst).x;
        self.data.iter_mut().for_each(|v| *v = *v * l);
        self.cat = Category::Illuminant;
        self
    }

    pub fn illuminance(&self, obs: &Observer) -> f64 {
        if self.cat == Category::Illuminant {
            (self.data * obs.data.column(1) * obs.lumconst).x
        } else {
            f64::NAN
        }
    }

    /// Downloads a spectrum
    pub async fn fetch(loc: &str) -> Result<Self, Box<dyn Error>> {
        let _url = Url::parse(loc)?;
        todo!()

    }

}

fn mixed_category(s1: &Spectrum, s2: &Spectrum) -> Category {
    if 
        (s1.cat == Category::Illuminant && s2.cat == Category::Colorant) ||
        (s1.cat == Category::Colorant && s2.cat == Category::Illuminant) ||
        (s1.cat == Category::Illuminant && s2.cat == Category::Filter) ||
        (s1.cat == Category::Filter && s2.cat == Category::Illuminant) {
        return Category::Stimulus
    } else if s1.cat == Category::Filter && s2.cat == Category::Filter {
            Category::Filter
    } else if s1.cat == Category::Colorant && s2.cat == Category::Colorant {
            Category::Colorant
    } else {
            Category::Unknown
    }
}

// Multiplication of Spectral, typically for cominations of an illuminant and a filter or colorant,
// or when combining multiple colorants or filters. Subtractive Mixing.
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
    ///     let mut led = Spectrum::led(550.0, 25.0);
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
            total: None
        }
    }
}

impl Mul<Spectrum> for f64 {
    /// Multiply a spectrum with a scalar f64 value.
    /// ```
    ///     use crate::colorimetry::Spectrum;
    ///     use approx::assert_ulps_eq;
    ///
    ///     let mut led = Spectrum::led(550.0, 25.0);
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
            total: None
        }
    }
}


// Addition of spectra, typically used for illuminant (multiple sources).
// Additive mixing
impl Add for Spectrum {
    type Output = Self;

    // multiply two cie spectra
    fn add(self, rhs: Self) -> Self::Output {
        let category = if self.cat == Category::Illuminant && rhs.cat == Category::Illuminant {
            Category::Illuminant
        } else {
            Category::Unknown
        };
        Self{
            cat: category,
            data: self.data.component_mul(&(rhs.data)),
            total: None
        }
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
    /// - for a colorant, this scales the total reflectivity.
    /// - for a filter, it changes its transmission.
    /// ```
    ///     use crate::colorimetry::Spectrum;
    ///     use approx::assert_ulps_eq;
    ///
    ///     let mut led = Spectrum::led(550.0, 25.0);
    ///     assert_ulps_eq!(led.irradiance(), 1.0, epsilon = 1E-10);
    ///
    ///     led *= 10.0;
    ///     assert_ulps_eq!(led.irradiance(), 10.0, epsilon = 1E-10);
    /// ```
    fn mul_assign(&mut self, rhs: f64) {
        self.data.iter_mut().for_each(|v| *v *= rhs);

    }
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
    use core::f64;

    use crate::spc::Spectrum;

    use crate::data::CIE1931;
    use approx::assert_ulps_eq;

    #[test]
    fn ee() {
        let [l, x, y ] = CIE1931.xyz(&&Spectrum::equal_energy().set_illuminance(&CIE1931, 100.0)).lxy();
        assert_ulps_eq!(l, 100.0, epsilon = f64::EPSILON);
        assert_ulps_eq!(x, 0.333_3, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.333_3, epsilon = 5E-5);
    }

    #[test]
    fn d65() {
        let [l, x, y ] = CIE1931.xyz(&Spectrum::d65().set_illuminance(&CIE1931, 100.0)).lxy();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(l, 100.0, epsilon = f64::EPSILON);
        assert_ulps_eq!(x, 0.312_72, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.329_03, epsilon = 5E-5);
    }

    #[test]
    fn d50() {
        let [l, x, y ] = CIE1931.xyz(&Spectrum::d50().set_illuminance(&CIE1931, 100.0)).lxy();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(l, 100.0, epsilon = 1E-8);
        assert_ulps_eq!(x, 0.345_67, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.358_51, epsilon = 5E-5);
    }

    #[test]
    fn a() {
        let [l, x, y ] = CIE1931.xyz(&Spectrum::a().set_illuminance(&CIE1931, 100.0)).lxy();
        // See table T3 CIE15:2004 (calculated with 5nm intervals, instead of 1nm, as used here)
        assert_ulps_eq!(l, 100.0, epsilon = f64::EPSILON);
        assert_ulps_eq!(x, 0.447_58, epsilon = 5E-5);
        assert_ulps_eq!(y, 0.407_45, epsilon = 5E-5);
    }
}