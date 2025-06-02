//! # Standard Observers
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

pub use observers::*;

use crate::{
    error::Error,
    illuminant::{CieIlluminant, Planck},
    lab::CieLab,
    math::LineAB,
    rgb::RgbSpace,
    spectrum::to_wavelength,
    spectrum::{Spectrum, NS, SPECTRUM_WAVELENGTH_RANGE},
    traits::{Filter, Light},
    xyz::XYZ,
};
use nalgebra::{Matrix3, SMatrix, Vector3};
use std::{ops::RangeInclusive, sync::OnceLock};
use strum_macros::EnumIter;
use wasm_bindgen::prelude::wasm_bindgen;

/**
   Light-weight identifier added to the `XYZ` and `RGB` datasets,
   representing the colorimetric standard observer used.

   No data included here, which would be the Rust way, but that does not work with wasm-bindgen.
   This can be directly used in JavaScript, and has the benefit to be just an index.
*/
#[cfg(not(feature = "supplemental-observers"))]
#[wasm_bindgen]
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug, EnumIter)]
pub enum Observer {
    #[default]
    Std1931,
}

#[cfg(feature = "supplemental-observers")]
#[wasm_bindgen]
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug, EnumIter)]
pub enum Observer {
    #[default]
    Std1931,
    Std1964,
    Std2015,
    Std2015_10,
}

impl Observer {
    /**
       Get a reference to the data for the specified `Observer`.
    */
    pub fn data(&self) -> &'static ObserverData {
        match self {
            Observer::Std1931 => &observers::CIE1931,
            #[cfg(feature = "supplemental-observers")]
            Observer::Std1964 => &observers::CIE1964,
            #[cfg(feature = "supplemental-observers")]
            Observer::Std2015 => &observers::CIE2015,
            #[cfg(feature = "supplemental-observers")]
            Observer::Std2015_10 => &observers::CIE2015_10,
        }
    }
}

/**
    A data structure to define Standard Observers, such as the CIE 1931 2º and
    the CIE 2015 standard observers.

    These are defined in the form of the three color matching functions,
    typically denoted by $\hat{x}(\lamda)$,$\hat{y}{\lambda}$, and $\hat{z}(\lambda)$.
    Traditionally, the CIE1931 Colorimetric Standard Observer is used almost exclusively,
    but is known to be not a good representation of human vision in the blue region of the
    spectrum. We also know now that the way you see color varies with age, and your healty,
    and that not everyone sees to same color.

    In this library colors are represented by spectral distributions, to allow color modelling
    with newer, and better standard observers, such as the CIE2015 Observer, derived from
    the sensitivities of the cones in the retina of your eye, the biological color receptors
    of light.

    It's main purpose is to calculate `XYZ` tristimulus values for a general stimulus,
    in from of a `Spectrum`.
*/
#[wasm_bindgen]
pub struct ObserverData {
    data: SMatrix<f64, 3, NS>,
    lumconst: f64,
    tag: Observer,
    d65: OnceLock<XYZ>,
    d50: OnceLock<XYZ>,

    /// The range of indices for which the spectral locus of this observer returns unique
    /// chromaticity coordinates. See documentation for the
    /// [`ObserverData::spectral_locus_wavelength_range`] method for details.
    spectral_locus_range: OnceLock<RangeInclusive<usize>>,
}

impl ObserverData {
    /// Creates a new `ObserverData` object, with the given color matching functions.
    ///
    /// Only visible to the crate itself since it cannot be used nicely from the outside
    /// (since the `tag` is not something anyone else can create new varians of).
    pub(crate) const fn new(tag: Observer, lumconst: f64, data: SMatrix<f64, 3, NS>) -> Self {
        Self {
            data,
            lumconst,
            tag,
            d65: OnceLock::new(),
            d50: OnceLock::new(),
            spectral_locus_range: OnceLock::new(),
        }
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
        let xyzn = light.xyzn(self.tag, None);
        if let Some(flt) = filter {
            let s = *light.spectrum() * *flt.spectrum();
            let xyz = self.xyz_from_spectrum(&s);
            let scale = 100.0 / xyzn.xyz.y;
            xyz * scale
        } else {
            xyzn
        }
    }

    /**
        Calculates Tristimulus valus, in form of an [XYZ] object of a general spectrum.
        If a reference white is given (rhs), it will copy its tristimulus value, and the spectrum
        is interpreted as a stimulus, being a combination of an illuminant with a colorant.
        If no reference white is given, the spectrum is interpreted as an illuminant.
        This method produces the raw XYZ data, not normalized to 100.0

    */
    pub fn xyz_from_spectrum(&self, spectrum: &Spectrum) -> XYZ {
        let xyz = self.data * spectrum.0 * self.lumconst;
        XYZ::from_vecs(xyz, self.tag)
    }

    /// Calculates the lumimous value or Y tristimulus value for a general spectrum.
    pub fn y_from_spectrum(&self, spectrum: &Spectrum) -> f64 {
        (self.data.row(1) * spectrum.0 * self.lumconst)[(0, 0)]
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
            .data
            .column(wavelength - SPECTRUM_WAVELENGTH_RANGE.start())
            .as_ref();
        Ok(XYZ::from_vecs(Vector3::new(x, y, z), self.tag))
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
    ///
    /// # Examples
    /// ```
    /// // Test data from [CIE 15:2018](https://cie.co.at/publications/colorimetry-4th-edition) for test data
    /// use colorimetry::{xyz::XYZ, illuminant::CieIlluminant};
    ///
    /// let cie1931_d65_xyz = colorimetry::observer::CIE1931.xyz_d65();
    /// approx::assert_ulps_eq!(cie1931_d65_xyz.values().as_ref(), [95.047, 100.0, 108.883].as_ref(), epsilon = 5E-2);
    /// ```
    pub fn xyz_d65(&self) -> XYZ {
        *self.d65.get_or_init(|| {
            self.xyz_from_spectrum(CieIlluminant::D65.illuminant().as_ref())
                .set_illuminance(100.0)
        })
    }

    /// XYZ tristimulus values for the CIE standard daylight illuminant D50 (buffered).
    ///
    /// # Examples
    /// ```ignore
    /// // Test data from [CIE 15:2018](https://cie.co.at/publications/colorimetry-4th-edition) for test data
    /// use colorimetry::{xyz::XYZ, illuminant::CieIlluminant};
    ///
    /// let cie1931_d50_xyz = colorimetry::observer::CIE1931.xyz_d50();
    /// approx::assert_ulps_eq!(cie1931_d50_xyz.values().as_ref(), [96.421, 100.0, 82.519].as_ref(), epsilon = 5E-2);
    ///
    /// # #[cfg(feature = "supplemental-observers")]
    /// # {
    /// let cie1964_d50_xyz = colorimetry::observer::CIE1964.xyz_d50();
    /// approx::assert_ulps_eq!(cie1964_d50_xyz.values().as_ref(), [96.720, 100.0, 81.427].as_ref(), epsilon = 5E-2);
    /// # }
    /// ```
    pub fn xyz_d50(&self) -> XYZ {
        *self.d50.get_or_init(|| {
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

    /// Calculates the L*a*b* CIELAB D65 values of a Colorant, using D65 as an illuminant.
    /// Accepts a Colorant Spectrum only.
    /// Returns f64::NAN's otherwise.
    pub fn lab_d65(&self, filter: &dyn Filter) -> CieLab {
        let xyz = self.xyz(&CieIlluminant::D65, Some(filter));
        let xyzn = self.xyz_d65();
        CieLab::from_xyz(xyz, xyzn).unwrap()
    }

    /// Calculates the L*a*b* CIELAB D50 values of a Colorant, using D65 as an illuminant.
    /// Accepts a Colorant Spectrum only.
    /// Returns f64::NAN's otherwise.
    pub fn lab_d50(&self, filter: &dyn Filter) -> CieLab {
        let xyz = self.xyz(&CieIlluminant::D50, Some(filter));
        let xyzn = self.xyz_d50();
        CieLab::from_xyz(xyz, xyzn).unwrap()
    }

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
    ///
    /// This range is computed on first access and then buffered for quick future access.
    pub fn spectral_locus_wavelength_range(&self) -> RangeInclusive<usize> {
        self.spectral_locus_range
            .get_or_init(|| {
                let min = self.spectral_locus_index_min() + *SPECTRUM_WAVELENGTH_RANGE.start();
                let max = self.spectral_locus_index_max() + *SPECTRUM_WAVELENGTH_RANGE.start();
                debug_assert!(min >= *SPECTRUM_WAVELENGTH_RANGE.start());
                debug_assert!(max <= *SPECTRUM_WAVELENGTH_RANGE.end());
                min..=max
            })
            .clone()
    }

    /// Calculates the RGB to XYZ matrix, for a particular color space.
    pub fn rgb2xyz(&self, rgbspace: &RgbSpace) -> Matrix3<f64> {
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
        self.rgb2xyz(&rgbspace).try_inverse().unwrap()
    }

    /// The index value of the blue spectral locus edge.
    ///
    /// Any further spectral locus points will hover around this edge, and will not have a unique wavelength.
    fn spectral_locus_index_min(&self) -> usize {
        const START: usize = 100;
        let spectral_locus_pos_start = self.spectral_locus_by_index(START).unwrap();
        let mut lp = LineAB::new(spectral_locus_pos_start, [0.33333, 0.33333]).unwrap();
        let mut m = START - 1;
        loop {
            let Some(spectral_locus_pos_m) = self.spectral_locus_by_index(m) else {
                break m + 1;
            };
            let l = LineAB::new(spectral_locus_pos_m, [0.33333, 0.33333]).unwrap();
            match (m, l.angle_diff(lp)) {
                (0, d) if d > -f64::EPSILON => break m + 1,
                (0, _) => break 0,
                (1.., d) if d > -f64::EPSILON => break m,
                _ => {
                    m -= 1;
                    lp = l;
                }
            }
        }
    }

    /// The index value of the red spectral locus edge.
    ///
    /// Any further spectral locus points will hover around this edge.
    fn spectral_locus_index_max(&self) -> usize {
        const START: usize = 300;
        let spectral_locus_pos_start = self.spectral_locus_by_index(START).unwrap();
        let mut lp = LineAB::new(spectral_locus_pos_start, [0.33333, 0.33333]).unwrap();
        let mut m = START + 1;
        loop {
            let Some(spectral_locus_pos_m) = self.spectral_locus_by_index(m) else {
                break m + 1;
            };
            let l = LineAB::new(spectral_locus_pos_m, [0.33333, 0.33333]).unwrap();
            match (m, l.angle_diff(lp)) {
                (400, d) if d < f64::EPSILON => break m - 1,
                (400, _) => break 400,
                (..400, d) if d < f64::EPSILON => break m - 1,
                _ => {
                    m += 1;
                    lp = l;
                }
            }
        }
    }

    /// Unrestricted, direct, access to the spectal locus data, in the form of
    /// chromaticity coordinates.
    fn spectral_locus_by_index(&self, i: usize) -> Option<[f64; 2]> {
        let &[x, y, z] = self.data.get((.., i))?.as_ref();
        let s = x + y + z;
        if s != 0.0 {
            Some([x / s, y / s])
        } else {
            None
        }
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
            .data
            .column_iter()
            .enumerate()
            .fold(Vector3::zeros(), |acc, (i, cmf)| {
                acc + cmf * f(i as f64 / (NS - 1) as f64)
            });
        XYZ {
            xyz,
            observer: self.tag,
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
    /// let rgb: [u8;3] = CIE1931.xyz_from_colorant_fn(&CieIlluminant::D65, |x|x).rgb(None).clamp().into();
    /// assert_eq!(rgb, [212, 171, 109]);
    /// ```
    /// Linear low pass filter, with a value of 1.0 for a wavelength of 380nm, and a value of 0.0 for 780nm,
    /// and converting the resulting value to RGB values.
    /// ```
    /// use colorimetry::prelude::*;
    /// let rgb: [u8;3] = CIE1931.xyz_from_colorant_fn(&CieIlluminant::D65, |x|1.0-x).rgb(None).clamp().into();
    /// assert_eq!(rgb, [158, 202, 237]);
    /// ```
    pub fn xyz_from_colorant_fn(&self, illuminant: &CieIlluminant, f: impl Fn(f64) -> f64) -> XYZ {
        let ill = illuminant.illuminant();
        let xyzn = self.xyz_cie_table(illuminant, None);
        let xyz = self
            .data
            .column_iter()
            .zip(ill.0 .0.iter())
            .enumerate()
            .fold(Vector3::zeros(), |acc, (i, (cmfi, sv))| {
                acc + cmfi * f(i as f64 / (NS - 1) as f64).clamp(0.0, 1.0) * *sv
            });
        let scale = 100.0 / xyzn.xyz.y;
        XYZ {
            xyz: xyz * self.lumconst * scale,
            observer: self.tag,
        }
    }
}

// JS-WASM Interface code
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl ObserverData {}

#[cfg(test)]
mod obs_test {

    use super::Observer;
    use crate::prelude::{CieIlluminant, CIE1931};
    use crate::rgb::RgbSpace;
    use crate::spectrum::SPECTRUM_WAVELENGTH_RANGE;
    use crate::xyz::XYZ;
    use approx::assert_ulps_eq;
    use strum::IntoEnumIterator as _;

    #[test]
    fn test_cie1931_spectral_locus_min_max() {
        let wavelength_lange = CIE1931.spectral_locus_wavelength_range();
        let chromaticity = CIE1931
            .xyz_at_wavelength(*wavelength_lange.start())
            .unwrap()
            .chromaticity();
        assert_ulps_eq!(chromaticity.x(), 0.17411, epsilon = 1E-5);
        assert_ulps_eq!(chromaticity.y(), 0.00496, epsilon = 1E-5);

        let chromaticity = CIE1931
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
            let wavelength_range = observer.data().spectral_locus_wavelength_range();

            // Basic sanity checking of the spectral locus wavelength range values
            assert!(wavelength_range.start() >= SPECTRUM_WAVELENGTH_RANGE.start());
            assert!(wavelength_range.end() <= SPECTRUM_WAVELENGTH_RANGE.end());
            assert!(wavelength_range.start() < wavelength_range.end());

            // Check that xyz_at_wavelength returns sane values in the allowed range
            for wavelength in wavelength_range {
                let xyz = observer.data().xyz_at_wavelength(wavelength).unwrap();
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
            eprintln!("Testing observer {:?}", observer);
            for wavelength in observer.data().spectral_locus_wavelength_range() {
                let xyz = observer.data().xyz_at_wavelength(wavelength).unwrap();
                for rgbspace in RgbSpace::iter() {
                    let _rgb = xyz.rgb(Some(rgbspace));
                }
            }
        }
    }

    #[test]
    // Data from http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
    // Lindbloom's values are reproduced with an accuracy of 3E-4, which is
    // a small, but significant difference.  This difference is caused by a difference in the display's white point,
    // due to wavelength domain differences.  Here we use a domain
    // from 380 to 780 with a step of 1 nanometer for the spectra, and in
    // specific for the color matching functions, as recommended by
    // CIE15:2004, and the color matching functions provided by the CIE in
    // their online dataset section. The spectra for the primaries are chosen to match the RGB primary values as given by the Color Space specifications.
    // The white point uses the standard's spectral distribution, as provided by the CIE, clipped to a domain from 380 to 780 nanometer.
    // See `colorimetry::data::D65`.
    fn test_rgb2xyz_cie1931() {
        let want = nalgebra::Matrix3::new(
            0.4124564, 0.3575761, 0.1804375, 0.2126729, 0.7151522, 0.0721750, 0.0193339, 0.1191920,
            0.9503041,
        );
        let got = CIE1931.rgb2xyz(&crate::rgb::RgbSpace::SRGB);
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
        let got = CIE1931.xyz2rgb(crate::rgb::RgbSpace::SRGB);
        approx::assert_ulps_eq!(want, got, epsilon = 3E-4);
    }

    /* not applicable anymore
    #[test]
    fn test_xyz_std_illuminants() {
        use crate::xyz::XYZ;
        use nalgebra as na;
        let XYZ {
            xyz,
            xyzn,
            observer,
        } = CIE1931.xyz_cie_table(&CieIlluminant::D65, Some(100.0));
        approx::assert_ulps_eq!(xyzn, na::Vector3::new(95.04, 100.0, 108.86), epsilon = 1E-2);
        assert!(xyz.is_none());
    }
     */

    #[test]
    fn test_planckian_locus() {
        // see https://www.waveformlighting.com/tech/calculate-cie-1931-xy-coordinates-from-cct
        // for test data (not clear what CMF domain they use)
        let xy = CIE1931
            .xyz_planckian_locus(3000.0)
            .chromaticity()
            .to_array();
        approx::assert_abs_diff_eq!(&xy.as_ref(), &[0.43693, 0.40407].as_ref(), epsilon = 2E-5);

        let xy = CIE1931
            .xyz_planckian_locus(6500.0)
            .chromaticity()
            .to_array();
        approx::assert_abs_diff_eq!(&xy.as_ref(), &[0.31352, 0.32363].as_ref(), epsilon = 6E-5);
    }

    #[test]
    fn test_xyz_from_illuminant_x_fn() {
        let xyz = CIE1931.xyz_from_colorant_fn(&CieIlluminant::D65, |_v| 1.0);
        let d65xyz = CIE1931.xyz_d65().xyz;
        approx::assert_ulps_eq!(
            xyz,
            crate::xyz::XYZ::from_vecs(d65xyz, crate::observer::Observer::Std1931)
        );
    }

    #[test]
    fn test_xyz_of_sample_with_standard_illuminant() {
        use crate::prelude::CieIlluminant::D65 as d65;
        let xyz = CIE1931
            .xyz(&d65, Some(&crate::colorant::Colorant::white()))
            .set_illuminance(100.0);
        approx::assert_ulps_eq!(xyz, CIE1931.xyz_from_colorant_fn(&d65, |_| 1.0));

        let xyz = CIE1931.xyz(&d65, Some(&crate::colorant::Colorant::black()));
        approx::assert_ulps_eq!(xyz, CIE1931.xyz_from_colorant_fn(&d65, |_| 0.0));
    }

    #[test]
    fn test_xyz_d65_d50() {
        let cie1931_d65_xyz = crate::observer::CIE1931.xyz_d65();
        approx::assert_ulps_eq!(
            cie1931_d65_xyz.values().as_ref(),
            [95.047, 100.0, 108.883].as_ref(),
            epsilon = 5E-2
        );

        let cie1931_d50_xyz = crate::observer::CIE1931.xyz_d50();
        approx::assert_ulps_eq!(
            cie1931_d50_xyz.values().as_ref(),
            [96.421, 100.0, 82.519].as_ref(),
            epsilon = 5E-2
        );
    }

    #[test]
    #[cfg(feature = "supplemental-observers")]
    fn test_xyz_d65_d50_cie1964() {
        let cie1964_d50_xyz = crate::observer::CIE1964.xyz_d50();
        approx::assert_ulps_eq!(
            cie1964_d50_xyz.values().as_ref(),
            [96.720, 100.0, 81.427].as_ref(),
            epsilon = 5E-2
        );

        let cie1964_d65_xyz = crate::observer::CIE1964.xyz_d65();
        approx::assert_ulps_eq!(
            cie1964_d65_xyz.values().as_ref(),
            [94.811, 100.0, 107.304].as_ref(),
            epsilon = 5E-2
        );
    }
}
