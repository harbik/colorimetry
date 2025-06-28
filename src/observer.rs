//! Standard Observers
//! ==================
//!
//! In color science, we model how humans perceive light using “standard observers”—mathematical averages of real people’s color responses. Here’s the gist:
//!
//! ## How We See Color
//! - Our eyes have three types of cones, each tuned to a different part of the spectrum:
//!   - **L-cones** (long-wavelength, “red” sensitive)  
//!   - **M-cones** (medium-wavelength, “green” sensitive)  
//!   - **S-cones** (short-wavelength, “blue” sensitive)  
//!
//! ## From Cone Sensitivities to XYZ
//!
//! Directly measuring each cone type is extremely challenging, so the CIE came up with a clever substitute:
//!
//! - **Color-matching experiments**  
//!   Observers tweak three narrow-band red, green, and blue lights until they visually match a test color.  
//! - **Negative matches**  
//!   Some test colors actually require “subtracting” a primary, which shows up as a negative match value.  
//! - **Imaginary primaries**  
//!   To eliminate negatives, the 1931 CIE team transformed those matches into three new “imaginary” primaries—mathematical constructs that guarantee all match values are zero or positive.  
//! - **Color-matching functions**  
//!   The non-negative weightings across wavelengths form the CIE 1931 curves x̅(λ), y̅(λ), and z̅(λ).  
//! - **XYZ tristimulus**  
//!   Finally, you integrate (dot-product) any light’s spectrum with these functions to get its standard X, Y, and Z values.
//!
//! ## Field of View Matters
//! - **2° Observer (CIE 1931)**  
//!   - Represents a small, foveal area of about 2°—think looking at a small color patch or point source.  
//!   - Based on experiments by Wright & Guild using narrow-band primaries and 17 observers.  
//! - **10° Observer (CIE 1964)**  
//!   - Covers a wider, more peripheral field (10°)—better for larger patches or immersive scenes.  
//!
//! ## Why New Observers?
//! - The original CIE 1931 functions work well across most of the spectrum but are less accurate in deep blue.  
//! - **CIE 2015 Observer** recalibrates those functions using updated cone-sensitivity data—especially improving blue-region accuracy.  
//! - Use the 2015 standard when you need the best possible match to human vision (e.g. advanced colorimetry, high-end display profiling).
//!
//! ## Which One to Use?
//! - **Small, detailed samples** → 2° (CIE 1931)  
//! - **Large fields or immersive scenes** → 10° (CIE 1964)  
//! - **Highest-accuracy applications** (especially blue-heavy content) → CIE 2015  
//!
//! These standard observers form the backbone of color spaces, chromaticity diagrams, and all standardized color measurements.  

mod observers;

mod optimal;
pub use optimal::OptimalColors;

mod spectral_locus;
pub use spectral_locus::SpectralLocus;

use crate::{
    cam::{CieCam02, CieCam16, ViewConditions},
    error::Error,
    illuminant::{CieIlluminant, Planck},
    lab::CieLab,
    rgb::RgbSpace,
    spectrum::{to_wavelength, Spectrum, NS, SPECTRUM_WAVELENGTH_RANGE},
    traits::{Filter, Light},
    xyz::{RelXYZ, XYZ},
};
use nalgebra::{Matrix3, SMatrix, Vector3};
use std::{fmt, ops::RangeInclusive, sync::OnceLock};

use strum::{AsRefStr, EnumIter};

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

/// Light-weight identifier added to the `XYZ` and `RGB` datasets,
///    representing the colorimetric standard observer used.
///
///    No data included here, which would be the Rust way, but that does not work with wasm-bindgen.
///    This can be directly used in JavaScript, and has the benefit to be just an index.
#[cfg(not(feature = "supplemental-observers"))]
#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, Debug, EnumIter, AsRefStr)]
#[non_exhaustive]
pub enum Observer {
    #[default]
    Cie1931,
}

#[cfg(feature = "supplemental-observers")]
#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, Debug, EnumIter, AsRefStr)]
#[non_exhaustive]
pub enum Observer {
    #[default]
    Cie1931,
    Cie1964,
    Cie2015,
    Cie2015_10,
}

impl Observer {
    /// Get a reference to the data for the specified `Observer`.
    fn data(&self) -> &'static ObserverData {
        match self {
            Observer::Cie1931 => &observers::CIE1931,
            #[cfg(feature = "supplemental-observers")]
            Observer::Cie1964 => &observers::CIE1964,
            #[cfg(feature = "supplemental-observers")]
            Observer::Cie2015 => &observers::CIE2015,
            #[cfg(feature = "supplemental-observers")]
            Observer::Cie2015_10 => &observers::CIE2015_10,
        }
    }

    /// Returns the name of the observer.
    pub fn name(&self) -> &'static str {
        self.data().name
    }

    /// Returns the spectral locus for this observer.
    ///
    /// The spectral locus is the boundary of the area of all physical colors in a chromiticity diagram.
    pub fn spectral_locus(&self) -> &SpectralLocus {
        self.data()
            .spectral_locus
            .get_or_init(|| SpectralLocus::new(*self))
    }

    pub fn rel_xyz(&self, light: &dyn Light, filter: &dyn Filter) -> RelXYZ {
        let xyzn = self.xyz_from_spectrum(&light.spectrum());
        let s = *light.spectrum() * *filter.spectrum();
        let xyz = self.xyz_from_spectrum(&s);
        let scale = 100.0 / xyzn.xyz.y;
        // unwrap OK as we are using only one observer (self) here
        RelXYZ::from_xyz(xyz * scale, xyzn * scale).unwrap()
    }

    /// Calulates Tristimulus values for an object implementing the [Light] trait, and an optional [Filter],
    /// filtering the light.
    ///
    /// The Light trait is implemented by [`CieIlluminant`] and [Illuminant](crate::illuminant::Illuminant).
    ///
    /// [`Colorant`](crate::colorant::Colorant) implments the [`Filter`] trait.
    /// [`Rgb`](crate::rgb::Rgb), which represents a display pixel, implements both in this library.
    /// As a light, it is the light emitted from the pixel, as a filter it is the RGB-composite
    /// filter which is applied to the underlying standard illuminant of color space.
    pub fn xyz(&self, light: &dyn Light, filter: Option<&dyn Filter>) -> XYZ {
        let xyzn = light.xyzn(self.data().tag, None);
        if let Some(flt) = filter {
            let s = *light.spectrum() * *flt.spectrum();
            let xyz = self.xyz_from_spectrum(&s);
            let scale = 100.0 / xyzn.xyz.y;
            xyz * scale
        } else {
            xyzn
        }
    }

    /// Calculates the L*a*b* CIELAB values for a light source and filter combination.
    /// This method is used to compute the color appearance of a light source
    /// when filtered by a colorant or filter.
    /// It returns the CIELAB values normalized to a white reference luminance of 100.0,
    ///
    /// # Arguments
    /// * `light` - A reference to an object implementing the [Light] trait, such as [`CieIlluminant`].
    /// * `filter` - A reference to an object implementing the [Filter] trait, such as `Colorant`.
    ///
    /// # Returns
    /// * `CieLab` - The computed CIELAB color representation for the light and filter combination.
    pub fn lab(&self, light: &dyn Light, filter: &dyn Filter) -> CieLab {
        let rxyz = self.rel_xyz(light, filter);
        // unwrap OK as we are using only one observer (self) here
        CieLab::from_xyz(rxyz)
    }

    /// Calculates the L*a*b* CIELAB D65 values of a Colorant, using D65 as an illuminant.
    /// Convenience method for `lab` with D65 illuminant, and using an illuminance value of 100.0.
    pub fn lab_d65(&self, filter: &dyn Filter) -> CieLab {
        self.lab(&crate::illuminant::D65, filter)
    }

    /// Calculates the L*a*b* CIELAB D50 values of a Colorant, using D50 as an illuminant.
    /// Convenience method for `lab` with D50 illuminant, and using an illuminance value of 100.0.
    pub fn lab_d50(&self, filter: &dyn Filter) -> CieLab {
        self.lab(&crate::illuminant::D50, filter)
    }

    /// Calculates the CIECAM16 color appearance model values for a light source and filter combination.
    /// This method is used to compute the color appearance of a light source
    /// when filtered by a colorant or filter, using the CIECAM16 model.
    /// It returns the CIECAM16 values normalized to a white reference luminance of 100.0.
    /// # Arguments
    /// * `light` - A reference to an object implementing the [Light] trait, such as [`CieIlluminant`].
    /// * `filter` - A reference to an object implementing the [Filter] trait, such as `Colorant`.
    /// * `vc` - The view conditions to use for the CIECAM16 calculation.
    /// # Returns
    /// * `CieCam16` - The computed CIECAM16 color appearance model representation for the light and filter combination.
    pub fn ciecam16(&self, light: &dyn Light, filter: &dyn Filter, vc: ViewConditions) -> CieCam16 {
        let rxyz = self.rel_xyz(light, filter);
        // unwrap OK as we are using only one observer (self) here
        CieCam16::from_xyz(rxyz, vc)
    }

    /// Calculates the CIECAM02 color appearance model values for a light source and filter combination.
    /// This method is used to compute the color appearance of a light source
    /// when filtered by a colorant or filter, using the CIECAM02 model.
    /// It returns the CIECAM02 values normalized to a white reference luminance of 100.0.
    /// # Arguments
    /// * `light` - A reference to an object implementing the [Light] trait, such as [`CieIlluminant`].         
    /// * `filter` - A reference to an object implementing the [Filter] trait, such as `Colorant`.
    /// * `vc` - The view conditions to use for the CIECAM02 calculation.
    /// # Returns
    /// * `CieCam02` - The computed CIECAM02 color appearance model representation for the light and filter combination.
    pub fn ciecam02(&self, light: &dyn Light, filter: &dyn Filter, vc: ViewConditions) -> CieCam02 {
        let rxyz = self.rel_xyz(light, filter);
        // unwrap OK as we are using only one observer (self) here
        CieCam02::from_xyz(rxyz, vc)
    }

    /// Calculates Tristimulus valus, in form of an [XYZ] object of a general spectrum.
    /// If a reference white is given (rhs), it will copy its tristimulus value, and the spectrum
    /// is interpreted as a stimulus, being a combination of an illuminant with a colorant.
    /// If no reference white is given, the spectrum is interpreted as an illuminant.
    /// This method produces the raw XYZ data, not normalized to 100.0
    pub fn xyz_from_spectrum(&self, spectrum: &Spectrum) -> XYZ {
        let xyz = self.data().data * spectrum.0 * self.data().lumconst;
        XYZ::from_vecs(xyz, self.data().tag)
    }

    /// Calculates the lumimous value or Y tristimulus value for a general spectrum.
    pub fn y_from_spectrum(&self, spectrum: &Spectrum) -> f64 {
        (self.data().data.row(1) * spectrum.0 * self.data().lumconst)[(0, 0)]
    }

    /// Returns the observer's color matching function (CMF) data as an [XYZ] tristimulus
    /// value for the given wavelength.
    ///
    /// This allows access to the underlying data for the entire range of wavelengths, 380-780nm.
    /// However, please read the documentation for the
    /// [`spectral_locus_wavelength_range`](Self::spectral_locus_wavelength_range) method on
    /// situations where you might not want to sample the full range.
    pub fn xyz_at_wavelength(&self, wavelength: usize) -> Result<XYZ, Error> {
        if !SPECTRUM_WAVELENGTH_RANGE.contains(&wavelength) {
            return Err(Error::WavelengthOutOfRange);
        };
        let &[x, y, z] = self
            .data()
            .data
            .column(wavelength - SPECTRUM_WAVELENGTH_RANGE.start())
            .as_ref();
        Ok(XYZ::from_vecs(Vector3::new(x, y, z), self.data().tag))
    }

    /// Calculates the relative XYZ tristimulus values of monochromatic stimuli.
    ///
    /// Monochromatic stimuli are pure spectral colors, like those seen when white light
    /// is dispersed by a prism. This method computes their XYZ values under a given illuminant.
    ///
    /// # Details
    /// - Computes XYZ values for wavelengths from 380nm to 780nm in 1nm steps
    /// - Multiplies observer color matching functions by illuminant spectrum
    /// - Normalizes results relative to the illuminant's white point
    /// - Scales output to 100 lux illuminance
    ///
    /// # Difference from Spectral Locus
    /// While the spectral locus shows only chromaticity coordinates (x,y) of pure spectral
    /// colors, this method provides full XYZ values including luminance information.
    ///
    /// # Implementation Notes
    /// - Each wavelength is treated as a monochromatic stimulus (delta function)
    /// - Results are typically low in magnitude due to the narrow bandwidth
    /// - Values are normalized relative to the illuminant's total energy
    ///
    /// # Parameters
    /// - `ref_white`: Reference illuminant (e.g., D65, D50) for normalization
    ///
    /// # Returns
    /// Vector of (wavelength, RelXYZ) pairs, where:
    /// - wavelength: nanometers (380-780nm)
    /// - RelXYZ: relative tristimulus values scaled to 100 lux
    pub fn monochromes(&self, ref_white: CieIlluminant) -> Vec<(usize, RelXYZ)> {
        let mut obs = self.data().data;
        let white = &ref_white.illuminant().as_ref().0;
        for r in 0..3 {
            for c in 0..NS {
                obs[(r, c)] *= white[c];
            }
        }
        let xyzn_vec = obs.column_sum();
        let xyzn = XYZ::from_vecs(xyzn_vec, self.data().tag);
        let mut v = Vec::with_capacity(NS);
        for w in SPECTRUM_WAVELENGTH_RANGE {
            let xyz = obs.column(w - SPECTRUM_WAVELENGTH_RANGE.start()).into();
            let rxyz = RelXYZ::from_vec(xyz, xyzn);
            v.push((w, rxyz.set_illuminance(100.0)));
        }
        v
    }

    pub fn trimmed_spectral_locus(&self, ref_white: CieIlluminant) -> Vec<(usize, RelXYZ)> {
        let sl_full = self.monochromes(ref_white);
        let valid_range = self.spectral_locus_wavelength_range();
        sl_full
            .into_iter()
            .filter(|(w, _)| valid_range.contains(w))
            .collect()
    }

    /// Tristimulus values for the Standard Illuminants in this library.
    ///
    /// Values are not normalized by default, unless an illuminance value is provided.
    ///
    /// TODO: buffer values
    pub fn xyz_cie_table(&self, std_illuminant: &CieIlluminant, illuminance: Option<f64>) -> XYZ {
        let xyz = self.xyz_from_spectrum(std_illuminant.illuminant().as_ref());
        if let Some(l) = illuminance {
            xyz.set_illuminance(l)
        } else {
            xyz
        }
    }

    /// XYZ tristimulus values for the CIE standard daylight illuminant D65 (buffered).
    pub fn xyz_d65(&self) -> XYZ {
        *self.data().d65.get_or_init(|| {
            self.xyz_from_spectrum(CieIlluminant::D65.illuminant().as_ref())
                .set_illuminance(100.0)
        })
    }

    /// XYZ tristimulus values for the CIE standard daylight illuminant D50 (buffered).
    pub fn xyz_d50(&self) -> XYZ {
        *self.data().d50.get_or_init(|| {
            self.xyz_from_spectrum(CieIlluminant::D50.illuminant().as_ref())
                .set_illuminance(100.0)
        })
    }

    /// Calculates XYZ tristimulus values for a Planckian emitter for this
    /// observer. The `to_wavelength`` function is used, as planck functions
    /// requires the wavelength to be in units of meters, and the
    /// `xyz_from_illuminant_as_fn` uses functions over a domain from 0.0 to
    /// 1.0.
    pub fn xyz_planckian_locus(&self, cct: f64) -> XYZ {
        let p = Planck::new(cct);
        self.xyz_from_fn(|l| p.at_wavelength(to_wavelength(l, 0.0, 1.0)))
    }

    /// The slope of the Plancking locus as a (dX/dT, dY/dT, dZ/dT) contained in
    /// a XYZ object.
    pub fn xyz_planckian_locus_slope(&self, cct: f64) -> XYZ {
        let p = Planck::new(cct);
        self.xyz_from_fn(|l| p.slope_at_wavelength(to_wavelength(l, 0.0, 1.0)))
    }

    /// Calculates the RGB to XYZ matrix, for a particular color space.
    pub fn rgb2xyz(&self, rgbspace: RgbSpace) -> Matrix3<f64> {
        let space = rgbspace.data();
        let mut rgb2xyz = Matrix3::from_iterator(
            space
                .primaries
                .iter()
                .flat_map(|s| self.xyz_from_spectrum(s).set_illuminance(1.0).values()),
        );
        let xyzw = self.xyz(&space.white, None).set_illuminance(1.0);
        let decomp = rgb2xyz.lu();
        // unwrap: only used with library color spaces
        let rgbw = decomp.solve(&xyzw.xyz).unwrap();
        for (i, mut col) in rgb2xyz.column_iter_mut().enumerate() {
            col *= rgbw[i];
        }
        rgb2xyz
    }

    /// Calculates the XYZ to RGB matrix, for a particular color space.
    pub fn xyz2rgb(&self, rgbspace: RgbSpace) -> Matrix3<f64> {
        // unwrap: only used with library color spaces
        self.rgb2xyz(rgbspace).try_inverse().unwrap()
    }

    //  pub fn rgb(&self, rgbspace: RgbSpace) -> Matrix3<f64> {
    //      self.rgb2xyz(rgbspace)
    //  }

    /// Returns the wavelength range (in nanometer) for the _horse shoe_,
    /// the boundary of the area of all physical colors in a chromiticity diagram,
    /// for this observer.
    ///
    /// Spectral locus points tend to freeze, or even fold back to lower wavelength
    /// values at the blue and red perimeter ends. This can be quite anoying, for
    /// example when trying to calculate dominant wavelength, or when creating
    /// plots.
    /// See Wikipedia's [CIE 1931 Color Space](https://en.wikipedia.org/wiki/CIE_1931_color_space).
    ///
    /// To help with the above problem, this method returns the wavelength range
    /// for which the spectral locus points are unique, meaning
    /// each wavelength has a chromaticity coordinate different from the wavelength
    /// below or above it.
    ///
    /// To get the tristimulus values of the spectral locus, use
    /// [`xyz_at_wavelength`](Self::xyz_at_wavelength).
    pub fn spectral_locus_wavelength_range(&self) -> RangeInclusive<usize> {
        self.data().spectral_locus_range.clone()
    }

    /// Calculates the XYZ tristimulus values for a spectrum defined by a function.
    ///
    /// The input function `f` should accept a floating-point value in the range `[0.0, 1.0]`,
    /// where `0.0` corresponds to a wavelength of 380 nm and `1.0` to 780 nm.
    /// The function will be called once for each wavelength step (401 times at 1 nm intervals).
    ///
    /// # Arguments
    /// * `f` - A function that takes a floating-point value in the range `[0.0, 1.0]` and returns
    ///   the spectral value at that wavelength, in units of watts per square meter per nanometer (W/m²/nm) for
    ///   illuminants, or Watts per square meter per steradian per nanometer (W/m²/sr/nm) for stimuli.
    ///
    /// # Notes
    /// - This method is used in the library to compute the Planckian locus (the color of blackbody
    ///   radiators), as described by Planck's law.
    /// - For colorants, use [`xyz_from_colorant_fn`](Self::xyz_from_colorant_fn).
    pub fn xyz_from_fn(&self, f: impl Fn(f64) -> f64) -> XYZ {
        let xyz = self
            .data()
            .data
            .column_iter()
            .enumerate()
            .fold(Vector3::zeros(), |acc, (i, cmf)| {
                acc + cmf * f(i as f64 / (NS - 1) as f64)
            });
        XYZ {
            xyz,
            observer: self.data().tag,
        }
    }

    /// Calculates XYZ tristimulus values for an analytic representation of a spectral distribution of
    /// a filter or a color patch, using a normalized wavelength domain ranging from a value of 0.0 to 1.0,
    /// illuminated with a standard illuminant.
    ///
    /// The spectral values should be defined within a range from 0.0 to 1.0, and are clamped otherwise.
    /// The resulting XYZ value will have relative Y values in the range from 0 to 100.0.
    ///
    /// # Examples
    /// A linear high pass filter, with a value of 0.0 for a wavelength of 380nm, and a value of 1.0 for 780nm,
    /// and converting the resulting value to RGB values.
    /// ```
    /// use colorimetry::prelude::*;
    /// use colorimetry::rgb::RgbSpace::SRGB;
    /// let rgb: [u8;3] = Cie1931.xyz_from_colorant_fn(&CieIlluminant::D65, |x|x).rgb(SRGB).clamp().into();
    /// assert_eq!(rgb, [212, 171, 109]);
    /// ```
    /// Linear low pass filter, with a value of 1.0 for a wavelength of 380nm, and a value of 0.0 for 780nm,
    /// and converting the resulting value to RGB values.
    /// ```
    /// use colorimetry::prelude::*;
    /// use colorimetry::rgb::RgbSpace::SRGB;
    /// let rgb: [u8;3] = Cie1931.xyz_from_colorant_fn(&CieIlluminant::D65, |x|1.0-x).rgb(SRGB).clamp().into();
    /// assert_eq!(rgb, [158, 202, 237]);
    /// ```
    pub fn xyz_from_colorant_fn(&self, illuminant: &CieIlluminant, f: impl Fn(f64) -> f64) -> XYZ {
        let ill = illuminant.illuminant();
        let xyzn = self.xyz_cie_table(illuminant, None);
        let xyz = self
            .data()
            .data
            .column_iter()
            .zip(ill.0 .0.iter())
            .enumerate()
            .fold(Vector3::zeros(), |acc, (i, (cmfi, sv))| {
                acc + cmfi * f(i as f64 / (NS - 1) as f64).clamp(0.0, 1.0) * *sv
            });
        let scale = 100.0 / xyzn.xyz.y;
        XYZ {
            xyz: xyz * self.data().lumconst * scale,
            observer: self.data().tag,
        }
    }
}

impl fmt::Display for Observer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.name().fmt(f)
    }
}

#[cfg(test)]
mod obs_test {

    #[cfg(feature = "supplemental-observers")]
    use super::Observer::Cie1964;
    use super::{Observer, Observer::Cie1931};
    use crate::{
        illuminant::CieIlluminant, rgb::RgbSpace, spectrum::SPECTRUM_WAVELENGTH_RANGE, xyz::XYZ,
    };
    use approx::assert_ulps_eq;
    use strum::IntoEnumIterator as _;

    #[test]
    fn test_rel_xyz() {
        use crate::colorant::Colorant;

        let obs = Cie1931;
        let light = CieIlluminant::D65;
        let filter = Colorant::gray(0.5);

        let rxyz = obs.rel_xyz(&light, &filter);
        assert_ulps_eq!(rxyz.xyz() * 2.0, rxyz.white_point(), epsilon = 1E-5);
        assert_ulps_eq!(rxyz.white_point(), obs.xyz_d65(), epsilon = 1E-5);
    }

    #[test]
    fn test_cie1931_spectral_locus_min_max() {
        let wavelength_lange = Cie1931.spectral_locus_wavelength_range();
        let chromaticity = Cie1931
            .xyz_at_wavelength(*wavelength_lange.start())
            .unwrap()
            .chromaticity();
        assert_ulps_eq!(chromaticity.x(), 0.17411, epsilon = 1E-5);
        assert_ulps_eq!(chromaticity.y(), 0.00496, epsilon = 1E-5);

        let chromaticity = Cie1931
            .xyz_at_wavelength(*wavelength_lange.end())
            .unwrap()
            .chromaticity();
        assert_ulps_eq!(chromaticity.x(), 0.73469, epsilon = 1E-5);
        assert_ulps_eq!(chromaticity.y(), 0.26531, epsilon = 1E-5);
    }

    /// Ensure all observers have sane spectral locus values
    #[test]
    fn test_spectral_locus_full() {
        for observer in Observer::iter() {
            let wavelength_range = observer.spectral_locus_wavelength_range();

            // Basic sanity checking of the spectral locus wavelength range values
            assert!(wavelength_range.start() >= SPECTRUM_WAVELENGTH_RANGE.start());
            assert!(wavelength_range.end() <= SPECTRUM_WAVELENGTH_RANGE.end());
            assert!(wavelength_range.start() < wavelength_range.end());

            // Check that xyz_at_wavelength returns sane values in the allowed range
            for wavelength in wavelength_range {
                let xyz = observer.xyz_at_wavelength(wavelength).unwrap();
                let chromaticity = xyz.chromaticity();
                assert!((0.0..=1.0).contains(&chromaticity.x()));
                assert!((0.0..=1.0).contains(&chromaticity.y()));

                // Check that all chromaticity coordinates on the spectral locus can be converted
                // to XYZ.
                let _xyz2 = XYZ::from_chromaticity(chromaticity, None, Some(observer)).unwrap();
            }
        }
    }

    /// Ensure XYZ values around all supported observer spectral locuses' can be converted to RGB
    #[test]
    fn test_spectral_locus_to_rgb() {
        for observer in Observer::iter() {
            eprintln!("Testing observer {observer:?}");
            for wavelength in observer.spectral_locus_wavelength_range() {
                let xyz = observer.xyz_at_wavelength(wavelength).unwrap();
                for rgbspace in RgbSpace::iter() {
                    let _rgb = xyz.rgb(rgbspace);
                }
            }
        }
    }

    #[test]
    // Data from http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html Lindbloom's values
    // are reproduced with an accuracy of 3E-4, which is a small, but significant difference.  This
    // difference is caused by a difference in the display's white point, due to wavelength domain
    // differences.  Here we use a domain from 380 to 780 with a step of 1 nanometer for the
    // spectra, and in specific for the color matching functions, as recommended by CIE15:2004, and
    // the color matching functions provided by the CIE in their online dataset section. The spectra
    // for the primaries are chosen to match the RGB primary values as given by the Color Space
    // specifications.  The white point uses the standard's spectral distribution, as provided by
    // the CIE, clipped to a domain from 380 to 780 nanometer.  See `colorimetry::data::D65`.
    fn test_rgb2xyz_cie1931() {
        let want = nalgebra::Matrix3::new(
            0.4124564, 0.3575761, 0.1804375, 0.2126729, 0.7151522, 0.0721750, 0.0193339, 0.1191920,
            0.9503041,
        );
        let got = Cie1931.rgb2xyz(crate::rgb::RgbSpace::SRGB);
        approx::assert_ulps_eq!(want, got, epsilon = 3E-4);
    }

    #[test]
    // Check the inverse transformation.
    // See comments at `test_rgb2xyz_cie1931`.
    fn test_xyz2rgb_cie1931() {
        let want = nalgebra::Matrix3::new(
            3.2404542, -1.5371385, -0.4985314, -0.9692660, 1.8760108, 0.0415560, 0.0556434,
            -0.2040259, 1.0572252,
        );
        let got = Cie1931.xyz2rgb(crate::rgb::RgbSpace::SRGB);
        approx::assert_ulps_eq!(want, got, epsilon = 3E-4);
    }

    #[test]
    fn test_planckian_locus() {
        // see https://www.waveformlighting.com/tech/calculate-cie-1931-xy-coordinates-from-cct
        // for test data (not clear what CMF domain they use)
        let xy = Cie1931
            .xyz_planckian_locus(3000.0)
            .chromaticity()
            .to_array();
        approx::assert_abs_diff_eq!(&xy.as_ref(), &[0.43693, 0.40407].as_ref(), epsilon = 2E-5);

        let xy = Cie1931
            .xyz_planckian_locus(6500.0)
            .chromaticity()
            .to_array();
        approx::assert_abs_diff_eq!(&xy.as_ref(), &[0.31352, 0.32363].as_ref(), epsilon = 6E-5);
    }

    #[test]
    fn test_xyz_from_illuminant_x_fn() {
        let xyz = Cie1931.xyz_from_colorant_fn(&CieIlluminant::D65, |_v| 1.0);
        let d65xyz = Cie1931.xyz_d65().xyz;
        approx::assert_ulps_eq!(
            xyz,
            crate::xyz::XYZ::from_vecs(d65xyz, crate::observer::Observer::Cie1931)
        );
    }

    #[test]
    fn test_xyz_of_sample_with_standard_illuminant() {
        use crate::prelude::CieIlluminant::D65 as d65;
        let xyz = Cie1931
            .xyz(&d65, Some(&crate::colorant::Colorant::white()))
            .set_illuminance(100.0);
        approx::assert_ulps_eq!(xyz, Cie1931.xyz_from_colorant_fn(&d65, |_| 1.0));

        let xyz = Cie1931.xyz(&d65, Some(&crate::colorant::Colorant::black()));
        approx::assert_ulps_eq!(xyz, Cie1931.xyz_from_colorant_fn(&d65, |_| 0.0));
    }

    #[test]
    fn test_xyz_d65_d50() {
        let cie1931_d65_xyz = Cie1931.xyz_d65();
        approx::assert_ulps_eq!(
            cie1931_d65_xyz.values().as_ref(),
            [95.047, 100.0, 108.883].as_ref(),
            epsilon = 5E-2
        );

        let cie1931_d50_xyz = Cie1931.xyz_d50();
        approx::assert_ulps_eq!(
            cie1931_d50_xyz.values().as_ref(),
            [96.421, 100.0, 82.519].as_ref(),
            epsilon = 5E-2
        );
    }

    #[test]
    #[cfg(feature = "supplemental-observers")]
    fn test_xyz_d65_d50_cie1964() {
        let cie1964_d50_xyz = Cie1964.xyz_d50();
        approx::assert_ulps_eq!(
            cie1964_d50_xyz.values().as_ref(),
            [96.720, 100.0, 81.427].as_ref(),
            epsilon = 5E-2
        );

        let cie1964_d65_xyz = Cie1964.xyz_d65();
        approx::assert_ulps_eq!(
            cie1964_d65_xyz.values().as_ref(),
            [94.811, 100.0, 107.304].as_ref(),
            epsilon = 5E-2
        );
    }

    #[test]
    fn test_spectral_locus() {
        use crate::illuminant::CieIlluminant;
        use crate::observer::Observer::Cie1931;

        let sl = Cie1931.monochromes(CieIlluminant::D65);
        // Check the first and last points of the spectral locus
        // Data obtained from spreadsheet using data directly downloaded from cie.co.at
        assert_eq!(sl.len(), 401);
        let xyz_first = sl[0].1.xyz().values();
        approx::assert_abs_diff_eq!(
            xyz_first.as_ref(),
            [6.46976E-04, 1.84445E-05, 3.05044E-03].as_ref(),
            epsilon = 1E-5
        );
        let xyz_last = sl[sl.len() - 1].1.xyz().values();
        approx::assert_abs_diff_eq!(
            xyz_last.as_ref(),
            [2.48982E-05, 8.99121E-06, 0.0].as_ref(),
            epsilon = 1E-5
        );
        // 550 nm
        let xyz_550 = sl[170].1.xyz().values();
        approx::assert_abs_diff_eq!(
            xyz_550.as_ref(),
            [0.4268018, 0.9796899, 0.0086158].as_ref(),
            epsilon = 1E-4
        );
    }

    #[test]
    fn test_as_ref_str() {
        // Ensure that the observer names are correct
        assert_eq!(Cie1931.as_ref(), "Cie1931");
        #[cfg(feature = "supplemental-observers")]
        {
            assert_eq!(Cie1964.as_ref(), "Cie1964");
            assert_eq!(Observer::Cie2015.as_ref(), "Cie2015");
            assert_eq!(Observer::Cie2015_10.as_ref(), "Cie2015_10");
        }
    }
}
