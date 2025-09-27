// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

//! Traits for spectrally-defined light sources and partially absorbing media
//! =========================================================================
//!
//! This module defines two key traits:
//!
//! - [`Light`]: a trait for anything that emits and can be described by a spectral power distribution.
//! - [`Filter`]: a trait for anything that (partially) absorbs light spectra, such as dyes and pigments.
//!
//! The [`Light`] trait is primarily used to represent illuminants or display pixel emissions,
//! providing tristimulus values for a given [`Observer`] through spectral integration.
//! This enables color calculations that are not limited to the CIE 1931 standard observer,
//! allowing modeling for animal vision, wide-field observers, or custom visual systems.
//!
//! You can implement `Light` for types that contain or can generate a [`Spectrum`].
//! For standard illuminants like D65, the trait is often implemented using lookup tables
//! or hard-coded spectral curves via the [`Illuminant`] type.
//!
//! ## Conversion
//! Any `Light` can be converted into an [`Illuminant`] using the `From` implementation.
//!
//! ## Related Types
//! - [`Spectrum`]: Represents sampled spectral data across wavelengths.
//! - [`Observer`]: Encapsulates color matching functions for various observers.
//! - [`Illuminant`]: Represents standard light sources such as D65 or A.
//!
//! [`Spectrum`]: crate::spectrum::Spectrum
//! [`Observer`]: crate::observer::Observer
//! [`Illuminant`]: crate::illuminant::Illuminant
use std::borrow::Cow;

use crate::{illuminant::Illuminant, observer::Observer, spectrum::Spectrum, xyz::XYZ};

/**
Spectral representation of Lights, typically in form of (standard) Illuminants.

Also allows to use lookup tristimulus values, such as the very common
[`D65`](crate::illuminant::CieIlluminant::D65) illuminant
(see [`CieIlluminant`](crate::illuminant::CieIlluminant) implemention).
Calculating them from a spectrum is the default implementation.
*/
pub trait Light {
    /// Calculates the tristimulus values of the light source, using the
    /// provided observer's color matching data, and an optional illuminance value.
    ///
    /// The illuminance value is optional, and if not provided, the actual luminous values in the
    /// spectrum are used.
    fn xyzn(&self, observer: Observer, y: Option<f64>) -> XYZ {
        let xyz = observer.xyz_from_spectrum(&self.spectrum());
        if let Some(illuminance) = y {
            xyz.set_illuminance(illuminance)
        } else {
            xyz
        }
    }

    /// Calculates the tristimulus values of the light source for a given observer,
    /// using an luminous value of 100.
    fn white_point(&self, observer: Observer) -> XYZ {
        self.xyzn(observer, Some(100.0))
    }

    #[cfg(feature = "cct")]
    fn cct(&self) -> Result<crate::illuminant::CCT, crate::Error> {
        let xyz = self.white_point(Observer::Cie1931);
        xyz.try_into()
    }

    fn spectrum(&self) -> Cow<'_, Spectrum>;
}

impl From<&dyn Light> for Illuminant {
    fn from(light: &dyn Light) -> Self {
        Illuminant(light.spectrum().into_owned())
    }
}

pub trait Filter {
    fn spectrum(&self) -> Cow<'_, Spectrum>;
}
