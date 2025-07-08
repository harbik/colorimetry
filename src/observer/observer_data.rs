use std::{ops::RangeInclusive, sync::OnceLock};

pub mod cie1931;
pub use cie1931::CIE1931;

pub mod cie1964;
pub use cie1964::CIE1964;

pub mod cie2015;
pub use cie2015::CIE2015;

pub mod cie2015_10;
pub use cie2015_10::CIE2015_10;

use nalgebra::SMatrix;

use crate::{
    observer::{Observer, SpectralLocus},
    spectrum::NS,
    xyz::XYZ,
};

///     A data structure to define Standard Observers, such as the CIE 1931 2º and
///     the CIE 2015 standard observers.
///
///     These are defined in the form of the three color matching functions,
///     typically denoted by $\hat{x}(\lamda)$,$\hat{y}{\lambda}$, and $\hat{z}(\lambda)$.
///     Traditionally, the CIE1931 Colorimetric Standard Observer is used almost exclusively,
///     but is known to be not a good representation of human vision in the blue region of the
///     spectrum. We also know now that the way you see color varies with age, and your healty,
///     and that not everyone sees to same color.
///
///     In this library colors are represented by spectral distributions, to allow color modelling
///     with newer, and better standard observers, such as the CIE2015 Observer, derived from
///     the sensitivities of the cones in the retina of your eye, the biological color receptors
///     of light.
///
///     It's main purpose is to calculate `XYZ` tristimulus values for a general stimulus,
///     in from of a `Spectrum`.
pub struct ObserverData {
    pub data: SMatrix<f64, 3, NS>,
    pub lumconst: f64,
    pub tag: Observer,
    pub name: &'static str,
    pub d65: OnceLock<XYZ>,
    pub d50: OnceLock<XYZ>,
    pub spectral_locus: OnceLock<SpectralLocus>,

    /// The range of indices for which the spectral locus of this observer returns unique
    /// chromaticity coordinates. See documentation for the
    /// [`ObserverData::spectral_locus_range`] method for details.
    pub spectral_locus_range: RangeInclusive<usize>,
}

impl ObserverData {
    /// Creates a new `ObserverData` object, with the given color matching functions.
    ///
    /// Only visible to the crate itself since it cannot be used nicely from the outside
    /// (since the `tag` is not something anyone else can create new varians of).
    const fn new(
        tag: Observer,
        name: &'static str,
        lumconst: f64,
        spectral_locus_range: RangeInclusive<usize>,
        data: SMatrix<f64, 3, NS>,
    ) -> Self {
        Self {
            data,
            lumconst,
            tag,
            name,
            d65: OnceLock::new(),
            d50: OnceLock::new(),
            spectral_locus_range,
            spectral_locus: OnceLock::new(),
        }
    }
}
