
use wasm_bindgen::prelude::wasm_bindgen;


/// The speed of light (m/s)
pub const C: f64 = 299792458.0;

/// Boltzmann constant (m<sup>2</sup> kg s<sup>-2</sup> K<sup>-1</sup>)
pub const KB: f64 = 1.3806485279E-23;

/// Planck constant (m<sup>2</sup> kg / s)
pub const H: f64 = 6.6260700408181E-34;

/// First radiation constant (W m<sup>2</sup>)
pub const C1: f64 = 2. * std::f64::consts::PI * H * C * C;

/// Second radiation constant (m K)
/// see <https://en.wikipedia.org/wiki/Planckian_locus#International_Temperature_Scale>
pub const C2: f64 = H * C / KB; // Now exact
pub const C2_NBS_1931: f64 = 1.435E-2; // A Illuminant
pub const C2_IPTS_1948: f64 = 1.4380E-2; // Illuminant series D
pub const C2_ITS_1968: f64 = 1.4388E-2;



/**
Planck with the second radiant constant as parameter.

This to be used for planckian radiators not in vacuum, or to calculate standard lights defined with
older values of this constant.

l: wavelength in meter
t: absolute temperature in Kelvin
c2: second radiative constant, in meter * Kelvin; can also be used to include refractive index, using c2 ::  c2 / n

*/
#[inline]
pub fn planck_c2(l: f64, t: f64, c2: f64) -> f64 {
    crate::physics::C1 / l.powi(5) / ((c2 / (l * t)).exp() - 1.0)
}

#[inline]
pub fn planck(l: f64, t: f64) -> f64 {
    planck_c2(l, t, crate::physics::C2)
}



/// Stefan-Boltzmann constant (W m<sup>-2</sup> K<sup>-4</sup>)
const SIGMA: f64 = 5.670_374_419_184E-8;

/// Stefan Boltzmann law: Blackbody's radiant emittance (W m<sup>-2</sup>), as function of its absolute
/// temperature (K).
#[inline]
#[wasm_bindgen(js_name= stefanBoltzmann)]
pub fn stefan_boltzmann(temperature: f64) -> f64 {
    SIGMA * temperature.powi(4)
}

const A: f64 = 1.11926158998; // scaling factor for power
const B: f64 = 1.08480681239; // scaling factor for width (l_w = B * l_fwhm)

/// LED Spectrum model
///
/// See Ohno, Spectral Design considerations for white LED Color Rendering, Optical Engineering 44(11), November 20005
/// Scale by spectralWidth
pub fn led_ohno(wl: f64, center: f64, width: f64) -> f64 {
    let width = B * width;
    let t = (wl - center) / (width);
    let g = libm::expm1(-(t.powi(2)))+1.0;
    (g + 2.0 * g.powi(5)) / (3.0 * A * width)
}

