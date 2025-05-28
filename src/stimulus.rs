use std::{
    borrow::Cow,
    iter::Sum,
    ops::{Deref, Mul},
};

use crate::{observer::ObserverData, rgb::Rgb, spectrum::Spectrum, traits::Light};

#[derive(Clone)]
pub struct Stimulus(pub(crate) Spectrum);

impl Deref for Stimulus {
    type Target = Spectrum;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Stimulus {
    /// Creates a new Stimulus from a spectrum.
    pub fn new(spectrum: Spectrum) -> Self {
        Self(spectrum)
    }

    /// Sets the luminance of the stimulus based on the observer data and a luminance value.
    pub fn set_luminance(mut self, obs: &ObserverData, luminance: f64) -> Self {
        let l = luminance / (obs.data.row(1) * self.0 .0 * obs.lumconst).x;
        self.0 .0.iter_mut().for_each(|v| *v *= l);
        self
    }

    /// A spectral composition of a display pixel, set to three sRGB color values.  The spectrum is
    /// a linear combination of the spectral primaries, which are Gaudssian filtered components in
    /// this library.
    pub fn from_srgb(r_u8: u8, g_u8: u8, b_u8: u8) -> Self {
        let rgb = Rgb::from_u8(
            r_u8,
            g_u8,
            b_u8,
            Some(crate::observer::Observer::Std1931),
            Some(crate::rgb::RgbSpace::SRGB),
        );
        rgb.into()
    }

    /// A spectral composition of a display pixel, set to three sRGB color values.  The spectrum is
    /// a linear combination of the spectral primaries, which are Gaudssian filtered components in
    /// this library.
    pub fn from_rgb(rgb: Rgb) -> Self {
        rgb.into()
    }
}

impl Light for Stimulus {
    fn spectrum(&self) -> Cow<Spectrum> {
        Cow::Borrowed(self)
    }
}

impl Sum for Stimulus {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut s = Spectrum::default();
        iter.for_each(|si| s += si.0);
        Stimulus(s)
    }
}

/// Spectral representation the color of a display pixel, described by a [`Rgb`]
/// instance.
///
/// It uses a linear combination of the spectral primaries as defined for a particular
/// [`RgbSpace``](crate::rgb::RgbSpace).
/// Most of the color spaces in this library use Daylight filtered Gaussian primaries,
/// but you can also use your own color space based on primaries measured by a spectrometer.
/// Spectral representations of pixels allow color matching for arbitrary observers,
/// not only the CIE 1931 standard observer.
impl From<Rgb> for Stimulus {
    fn from(rgb: Rgb) -> Self {
        let prim = &rgb.space.data().primaries;
        let rgb2xyz = rgb.observer.data().rgb2xyz(&rgb.space);
        let yrgb = rgb2xyz.row(1);
        rgb.rgb
            .iter()
            .zip(yrgb.iter())
            .zip(prim.iter())
            .map(|((v, w), s)| *v * *w * s.clone())
            .sum()
    }
}

impl Mul<f64> for Stimulus {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self(self.0 * rhs)
    }
}

impl Mul<Stimulus> for f64 {
    type Output = Stimulus;

    fn mul(self, rhs: Stimulus) -> Self::Output {
        Stimulus(self * rhs.0)
    }
}
