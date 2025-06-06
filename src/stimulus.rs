//! Spectral Power Distributions for Visible Light Stimuli
//! ======================================================
//!
//! This module defines the [`Stimulus`] struct, which represents the spectral power
//! distribution of a visible light stimulus—such as a pixel on a screen—using a [`Spectrum`].
//!
//! A `Stimulus` is typically constructed from RGB values using `from_srgb` or `from_rgb`,
//! mapping them to a linear combination of spectral primaries defined by a color space
//! (e.g., sRGB with Gaussian-filtered components). These spectral primaries allow for
//! physically-informed color calculations that are observer-aware.
//!
//! Unlike traditional RGB or XYZ values, a `Stimulus` retains full spectral detail,
//! enabling colorimetric computations for arbitrary [`Observer`] types,
//! not just the standard CIE 1931 observer.
//!
//! ## Features
//! - Convert RGB values into spectral representations.
//! - Scale to a target luminance using real observer sensitivity curves.
//! - Integrates with the [`Light`] trait for uniform handling of illuminants and reflectances.
//!
//! ## When to Use `Stimulus`
//! Use `Stimulus` when:
//! - You need to model light with spectral fidelity.
//! - You want to convert RGB colors to spectra for metamerism studies or non-standard observers.
//! - You need to simulate how different humans perceive color.
//!
//! [`Spectrum`]: crate::spectrum::Spectrum
//! [`Observer`]: crate::observer::Observer
//! [`Light`]: crate::traits::Light

use std::{
    borrow::Cow,
    iter::Sum,
    ops::{Deref, Mul},
};

use crate::{observer::Observer, rgb::Rgb, spectrum::Spectrum, traits::Light};

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
    pub fn set_luminance(mut self, obs: Observer, luminance: f64) -> Self {
        let y = obs.y_from_spectrum(self.as_ref());
        let l = luminance / y;
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
            Some(crate::observer::Observer::Cie1931),
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
        let rgb2xyz = rgb.observer.rgb2xyz(rgb.space);
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

impl AsRef<Spectrum> for Stimulus {
    fn as_ref(&self) -> &Spectrum {
        &self.0
    }
}
