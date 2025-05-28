//! # Planck Radiation model and Physical Constants
//!

/// The speed of light (m/s)
const C: f64 = 299792458.0;

/// Boltzmann constant (m<sup>2</sup> kg s<sup>-2</sup> K<sup>-1</sup>)
const KB: f64 = 1.3806485279E-23;

/// Planck constant (m<sup>2</sup> kg / s)
const H: f64 = 6.6260700408181E-34;

/// First radiation constant (W m<sup>2</sup>)
const C1: f64 = 2. * std::f64::consts::PI * H * C * C;

/// Variants of the second radiation constant \(c_2\) used in Planck’s law,
/// reflecting historical definitions and the current exact definition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SecondRadiationConstant {
    /// Current exact definition: \(c_2 = \frac{h \, c}{k_B}\)
    Exact,
    /// Value used in the definition of the A Illuminant (NBS 1931)
    Nbs1931,
    /// Value used in the D Illuminant series (IPTS 1948)
    Ipts1948,
    /// Value used in the D Illuminant series (ITS 1968)
    Its1968,
}

impl SecondRadiationConstant {
    /// Returns the numerical value of the chosen constant (in m·K).
    pub fn value(self) -> f64 {
        match self {
            SecondRadiationConstant::Exact => H * C / KB,
            SecondRadiationConstant::Nbs1931 => 1.435e-2,
            SecondRadiationConstant::Ipts1948 => 1.4380e-2,
            SecondRadiationConstant::Its1968 => 1.4388e-2,
        }
    }
}

/// Planck's law describes the spectral radiance of a black body at a given temperature.  This
/// struct encapsulates the temperature of the black body, with unit of Kelvin, and provides methods
/// to calculate
pub struct Planck(f64);

impl Planck {
    /// Create a new Planck instance with the given temperature with a unit of Kelvin.
    pub fn new(temperature: f64) -> Self {
        Planck(temperature)
    }

    /// Calculate the spectral radiance at a given wavelength.
    /// This is based on Planck's law, which describes the spectral radiance of a black body at a given temperature.
    /// # Arguments
    /// * `wavelength` - The wavelength in meters at which to calculate the spectral radiance.
    /// # Returns
    /// The spectral radiance in watts per square meter per steradian per meter (W·m<sup>-2</sup>·sr<sup>-1</sup>·m<sup>-1</sup>).
    pub fn at_wavelength(&self, wavelength: f64) -> f64 {
        let t = self.0;
        let c2 = SecondRadiationConstant::Exact.value();
        C1 / wavelength.powi(5) / ((c2 / (wavelength * t)).exp() - 1.0)
    }

    /// Calculate the slope of the spectral radiance with respect to temperature.
    /// This is the first derivative of Planck's law.
    /// # Arguments
    /// * `wavelength` - The wavelength in meters at which to calculate the slope.
    /// # Returns
    pub fn slope_at_wavelength(&self, wavelength: f64) -> f64 {
        let t = self.0;
        let c2 = SecondRadiationConstant::Exact.value();
        let c3 = C1 * c2 / t.powi(2);
        let e = (c2 / (wavelength * t)).exp();
        c3 / wavelength.powi(6) * e / (e - 1.0).powi(2)
    }

    pub fn planck_with_legacy_c2(&self, wavelength: f64, c2: SecondRadiationConstant) -> f64 {
        let t = self.0;
        let c2_value = c2.value();
        C1 / wavelength.powi(5) / ((c2_value / (wavelength * t)).exp() - 1.0)
    }

    pub fn slope_with_legacy_c2(&self, wavelength: f64, c2: SecondRadiationConstant) -> f64 {
        let t = self.0;
        let c2_value = c2.value();
        let c3 = C1 * c2_value / t.powi(2);
        let e = (c2_value / (wavelength * t)).exp();
        c3 / wavelength.powi(6) * e / (e - 1.0).powi(2)
    }

    pub fn curvature_with_legacy_c2(&self, wavelength: f64, c2: SecondRadiationConstant) -> f64 {
        let t = self.0;
        let c2_value = c2.value();
        let e = (c2_value / (wavelength * t)).exp();
        self.slope_with_legacy_c2(wavelength, c2) / t
            * (c2_value / (wavelength * t) * (e + 1.0) / (e - 1.0) - 2.0)
    }

    /// Calculate the total radiant emittance of a black body at the given temperature.
    /// This is based on the Stefan-Boltzmann law, which states that the total energy radiated per unit surface area
    /// of a black body is proportional to the fourth power of its absolute temperature.
    pub fn stefan_boltzmann(&self) -> f64 {
        let temperature = self.0;
        const SIGMA: f64 = 5.670374419184E-8; // Stefan-Boltzmann constant (W m<sup>-2</sup> K<sup>-4</sup>)
        SIGMA * temperature.powi(4)
    }
}
