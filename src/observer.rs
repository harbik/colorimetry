/*!

# Standard Observers

- Human color perception
    - begins in the retina, triggered by light stimuli detected by the **retina’s L, M, and S cones**,
      each sensitive to different wavelength ranges.
      - The **L cone** is most sensitive to long wavelengths (red),
      - the **M cone** to medium wavelengths (green),
      - and the **S cone** to short wavelengths (blue).

- Visual responsivity varies with the angular extent of the field of view.
  - The **CIE 1931 Standard Observer** defines color matching functions for a **2-degree** field of
    view, representing central (foveal) vision.  This is typically used for viewing **small objects**,
    such as point light sources or fine details.
  - The **CIE 1964 Standard Observer** extends the field of view to **10 degrees**, capturing a
    broader area of the retina, including more peripheral cones.  It is recommended for **larger
    objects** or extended fields, such as color patches, or when measuring color in more immersive
    viewing conditions.
- Cone response differ across individuals, leading to variations in color perception.
  - **Genetic variation**, especially on the **X chromosome**, causes differences in cone types and
    densities—explaining individual variations and color vision deficiencies (e.g., red-green color
    blindness).
  - Additional variation comes from the **optical media of the eye** (lens, cornea, macula), which
    filter incoming light and change with age and health.
- Cone responses form the physiological basis of **tristimulus values (X, Y, Z)**, first defined by the **CIE**
  in 1931.
  - Color matching experiments using **bipartite fields** and narrow-band **RGB primaries** were used
    to measure how observers match colors across the visible spectrum.
  - Data from **Wright** and **Guild**, using a total of 17 male observers, were combined to define
    average **RGB color matching functions**.
  - To eliminate negative values, the RGB functions were transformed into the **XYZ space** using
    imaginary primaries.
  - The Y function was aligned with the **CIE 1924 photopic luminosity function**.
  - The resulting **color matching functions**—𝑥̅(λ), 𝑦̅(λ), and 𝑧̅(λ)— define the standard observer’s response.
  - These functions are the foundation of **chromaticity diagrams** like the **CIE 1931 (x, y)**
    space, used to visualize color perception across the visible spectrum.
- The **CIE 1931 Standard Observer** is the most widely used, but it has limitations.
  - It does not accurately represent color perception in the blue region of the spectrum.
  - The **CIE 2015 Standard Observer** was developed to address these limitations, providing a more
    accurate representation of human color vision, especially in the blue region.
  - It is based on the **CIE 2015 color matching functions**, which are derived from the spectral
    sensitivities of the cones in the retina.
  - The CIE 2015 Standard Observer is recommended for applications requiring high color accuracy,
    such as colorimetry, color science, and display technology.
  - It is also used in the **CIE 2015 color space**, which is a more accurate representation of
    color perception than the CIE 1931 color space.

 */

use crate::{
    colorant::Colorant,
    error::CmtError,
    geometry::LineAB,
    lab::CieLab,
    physics::{planck, planck_slope, to_wavelength},
    rgbspace::RgbSpace,
    spectrum::{Spectrum, NS, SPECTRUM_WAVELENGTH_RANGE},
    std_illuminants::StdIlluminant,
    traits::{Filter, Light},
    xyz::XYZ,
};
use nalgebra::{Matrix3, SMatrix, Vector3};
use std::{
    borrow::{Borrow, Cow},
    sync::OnceLock,
};
use wasm_bindgen::{convert::IntoWasmAbi, prelude::wasm_bindgen};

/**
   Light-weight identifier added to the `XYZ` and `RGB` datasets,
   representing the colorimetric standard observer used.

   No data included here, which would be the Rust way, but that does not work with wasm-bindgen.
   This can be directly used in JavaScript, and has the benefit to be just an index.
*/
#[cfg(not(feature = "supplemental-observers"))]
#[wasm_bindgen]
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub enum Observer {
    #[default]
    Std1931,
}

#[cfg(feature = "supplemental-observers")]
#[wasm_bindgen]
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
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
            Observer::Std1931 => &crate::data::observers::CIE1931,
            #[cfg(feature = "supplemental-observers")]
            Observer::Std1964 => &crate::data::observers::CIE1964,
            #[cfg(feature = "supplemental-observers")]
            Observer::Std2015 => &crate::data::observers::CIE2015,
            #[cfg(feature = "supplemental-observers")]
            Observer::Std2015_10 => &crate::data::observers::CIE2015_10,
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
    pub(crate) data: SMatrix<f64, 3, NS>,
    pub(crate) lumconst: f64,
    pub(crate) tag: Observer,
}

impl ObserverData {
    /// Calulates Tristimulus values for an object implementing the [Light] trait, and an optional [Filter],
    /// filtering the light.
    ///
    /// The Light trait is implemented by [`StdIlluminant`] and [Illuminant](crate::illuminant::Illuminant).
    ///
    /// [`Colorant`] implments the [`Filter`] trait.
    /// [`RGB`](crate::rgb::RGB), which represents a display pixel, implements both in this library.
    /// As a light, it is the light emitted from the pixel, as a filter it is the RGB-composite
    /// filter which is applied to the underlying standard illuminant of color space.
    pub fn xyz(&self, light: &dyn Light, filter: Option<&dyn Filter>) -> XYZ {
        // let xyzn = self.xyz_cie_table(illuminant, None);
        let xyzn = light.xyzn(self.tag, None);
        let xyz = if let Some(flt) = filter {
            let s = *light.spectrum() * *flt.spectrum();
            self.xyz_from_spectrum(&s, Some(xyzn))
        } else {
            xyzn
        };
        xyz.set_illuminance(100.0)
    }

    /**
        Calculates Tristimulus valus, in form of an [XYZ] object of a general spectrum.
        If a reference white is given (rhs), it will copy its tristimulus value, and the spectrum
        is interpreted as a stimulus, being a combination of an illuminant with a colorant.
        If no reference white is given, the spectrum is interpreted as an illuminant.
        This method produces the raw XYZ data, not normalized to 100.0

    */
    pub fn xyz_from_spectrum(&self, spectrum: &Spectrum, rhs: Option<XYZ>) -> XYZ {
        let xyz = self.data * spectrum.0 * self.lumconst;
        if let Some(xyz0) = rhs {
            // A illuminant/colorant
            XYZ::from_vecs(xyz0.xyzn, Some(xyz), self.tag)
        } else {
            // illuminant only
            XYZ::from_vecs(xyz, None, self.tag)
        }
    }

    /**
        Tristimulus Values for the Standard Illuminants in this library.

        Values are calculated on first use, and are not normalized by default, unless an illuminous
        value is provided, in case they are.
    */
    pub fn xyz_cie_table(&self, std_illuminant: &StdIlluminant, illuminance: Option<f64>) -> XYZ {
        const XYZ_STD_ILLUMINANTS_LEN: usize = 64;
        static XYZ_STD_ILLUMINANTS: [OnceLock<XYZ>; XYZ_STD_ILLUMINANTS_LEN] =
            [const { OnceLock::new() }; XYZ_STD_ILLUMINANTS_LEN];

        let xyz = *XYZ_STD_ILLUMINANTS[*std_illuminant as usize]
            .get_or_init(|| self.xyz_from_spectrum(std_illuminant.illuminant(), None));
        if let Some(l) = illuminance {
            xyz.set_illuminance(l)
        } else {
            xyz
        }
    }

    /// XYZ tristimulus values for the CIE standard daylight illuminant D65.
    /// The values are calculated on first use.
    pub fn xyz_d65(&self) -> XYZ {
        self.xyz_cie_table(&StdIlluminant::D65, Some(100.0))
    }

    /// XYZ tristimulus values for the CIE standard daylight illuminant D50.
    /// The values are calculated on first use.
    pub fn xyz_d50(&self) -> XYZ {
        self.xyz_cie_table(&StdIlluminant::D50, Some(100.0))
    }

    /**
        Calculates XYZ tristimulus values for an analytic representation of a spectral distribution of
        a filter or a color patch, using a normalized wavelength domain ranging from a value of 0.0 to 1.0,
        illuminated with a standard illuminant.

        The spectral values should be defined within a range from 0.0 to 1.0, and are clamped otherwise.
        The resulting XYZ value will have relative Y values in the range from 0 to 100.0,
        and yw is set to a value of 100.0.

        # Examples
        Linear high pass filter, with a value of 0.0 for a wavelength of 380nm, and a value of 1.0 for 780nm,
        and converting the resulting value to RGB values.
        ```
            use colorimetry::prelude::*;
            let rgb: [u8;3] = CIE1931.xyz_from_std_illuminant_x_fn(&StdIlluminant::D65, |x|x).rgb(None).into();
            assert_eq!(rgb, [212, 171, 109]);
        ```
        Linear low pass filter, with a value of 1.0 for a wavelength of 380nm, and a value of 0.0 for 780nm,
        and converting the resulting value to RGB values.
        ```
            use colorimetry::prelude::*;
            let rgb: [u8;3] = CIE1931.xyz_from_std_illuminant_x_fn(&StdIlluminant::D65, |x|1.0-x).rgb(None).into();
            assert_eq!(rgb, [158, 202, 237]);
        ```

    */
    pub fn xyz_from_std_illuminant_x_fn(
        &self,
        illuminant: &StdIlluminant,
        f: impl Fn(f64) -> f64,
    ) -> XYZ {
        let s = illuminant.illuminant();
        let xyzn = self.xyz_cie_table(illuminant, None).xyzn;
        let xyz =
            self.data.column_iter().zip(s.0 .0.iter()).enumerate().fold(
                Vector3::zeros(),
                |acc, (i, (cmfi, sv))| {
                    acc + cmfi * f(i as f64 / (NS - 1) as f64).clamp(0.0, 1.0) * *sv
                },
            );
        XYZ {
            xyz: Some(xyz * self.lumconst),
            xyzn,
            observer: self.tag,
        }
        .set_illuminance(100.0)
    }

    /**
        Calculates XYZ tristimulus values for an illuminant with its spectral distribution
        described by a function, defined over a domain from 0.0 to 1.0, with 0.0 corresponding to
        a wavelength of 380nm, and 1.0 to a wavelength of 780nm.

        It is mainly used in this library to calculate the Planckian locus, which is described by
        Planck's law.  The resulting XYZ value will be normalized to hava a Y value of 100.0
        and yw is set to None.
    */
    pub fn xyz_fn_illuminant(&self, f: impl Fn(f64) -> f64) -> XYZ {
        let xyzn = self
            .data
            .column_iter()
            .enumerate()
            .fold(Vector3::zeros(), |acc, (i, cmf)| {
                acc + cmf * f(i as f64 / (NS - 1) as f64)
            });
        XYZ {
            xyz: None,
            xyzn,
            observer: self.tag,
        }
    }

    /// Calculates XYZ tristimulus values for a Planckian emitter for this
    /// observer. The `to_wavelength`` function is used, as planck functions
    /// requires the wavelength to be in units of meters, and the
    /// `xyz_from_illuminant_as_fn` uses functions over a domain from 0.0 to
    /// 1.0.
    pub fn xyz_planckian_locus(&self, cct: f64) -> XYZ {
        self.xyz_fn_illuminant(|l| planck(to_wavelength(l, 0.0, 1.0), cct))
    }

    /// The slope of the Plancking locus as a (dX/dT, dY/dT, dZ/dT) contained in
    /// a XYZ object.
    pub fn xyz_planckian_locus_slope(&self, cct: f64) -> XYZ {
        self.xyz_fn_illuminant(|l| planck_slope(to_wavelength(l, 0.0, 1.0), cct))
    }

    /// Calculates the L*a*b* CIELAB D65 values of a Colorant, using D65 as an illuminant.
    /// Accepts a Colorant Spectrum only.
    /// Returns f64::NAN's otherwise.
    pub fn lab_d65(&self, filter: &dyn Filter) -> CieLab {
        let xyz0 = self.xyz(&StdIlluminant::D65, Some(filter));
        xyz0.try_into().unwrap()
    }

    /// Calculates the L*a*b* CIELAB D50 values of a Colorant, using D65 as an illuminant.
    /// Accepts a Colorant Spectrum only.
    /// Returns f64::NAN's otherwise.
    pub fn lab_d50(&self, filter: &dyn Filter) -> CieLab {
        let xyz0 = self.xyz(&StdIlluminant::D50, Some(filter));
        xyz0.try_into().unwrap()
    }

    /// Calculate the Spectral Locus, or (x,y) coordinates of the _horse shoe_,
    /// the boundary of area of all physical colors in a chromiticity diagram,
    /// as XYZ tristimulus values.  Tristimulus values are returned here,
    /// instead of chromaticity xy coordinates, to make use of the available XYZ
    /// transforms.
    ///
    /// This function limits the input values to produce unique chromaticity
    /// values only.  Spectral locus points tend to freeze, or even fold back to
    /// lower wavelength values at the blue and red perimeter ends.  This can be
    /// quite anoying, for example when trying to calculate dominant wavelength,
    /// or when creating plots.  To get the allowed range, use the
    /// `spectral_wavelength_min` and `spectral_wavelength_max` methods.
    ///
    /// See Wikipedia's [CIE 1931 Color Space](https://en.wikipedia.org/wiki/CIE_1931_color_space).
    pub fn spectral_locus_by_nm(&self, l: usize) -> Result<XYZ, CmtError> {
        let min = self.spectral_locus_nm_min();
        let max = self.spectral_locus_nm_max();
        if !SPECTRUM_WAVELENGTH_RANGE.contains(&l) {
            return Err(CmtError::WavelengthOutOfRange);
        };
        if l < min || l > max {
            Err(CmtError::NoUniqueSpectralLocus(min, max))
        } else {
            let &[x, y, z] = self.data.column(l - 380).as_ref();
            Ok(XYZ::from_vecs(Vector3::new(x, y, z), None, self.tag))
        }
    }

    /// Unrestricted, direct, access to the spectal locus data.
    /// To get unique values only please use the `spectral_locus_by_nm` function.
    pub fn spectral_locus_by_index(&self, i: usize) -> [f64; 2] {
        let &[x, y, z] = self.data.column(i).as_ref();
        let s = x + y + z;
        [x / s, y / s]
    }

    /// The index value of the blue spectral locus edge.
    ///
    /// Any further spectral locus points will hover around this edge, and will not have a unique wavelength.
    pub fn spectral_locus_index_min(&self) -> usize {
        static MIN: OnceLock<usize> = OnceLock::new();
        *MIN.get_or_init(|| {
            const START: usize = 100;
            let mut lp =
                LineAB::new(self.spectral_locus_by_index(START), [0.33333, 0.33333]).unwrap();
            let mut m = START - 1;
            loop {
                let l = LineAB::new(self.spectral_locus_by_index(m), [0.33333, 0.33333]).unwrap();
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
        })
    }

    pub fn spectral_locus_nm_min(&self) -> usize {
        self.spectral_locus_index_min() + 380
    }

    /// The index value of the red spectral locus edge.
    ///
    /// Any further spectral locus points will hover around this edge.
    pub fn spectral_locus_index_max(&self) -> usize {
        static MAX: OnceLock<usize> = OnceLock::new();
        *MAX.get_or_init(|| {
            const START: usize = 300;
            let mut lp =
                LineAB::new(self.spectral_locus_by_index(START), [0.33333, 0.33333]).unwrap();
            let mut m = START + 1;
            loop {
                let l = LineAB::new(self.spectral_locus_by_index(m), [0.33333, 0.33333]).unwrap();
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
        })
    }

    pub fn spectral_locus_nm_max(&self) -> usize {
        self.spectral_locus_index_max() + 380
    }

    /// Calculates the RGB to XYZ matrix, for a particular color space.
    pub fn rgb2xyz(&self, rgbspace: &RgbSpace) -> Matrix3<f64> {
        let space = rgbspace.data();
        let mut rgb2xyz = Matrix3::from_iterator(space.primaries.iter().flat_map(|s| {
            self.xyz_from_spectrum(s, None)
                .set_illuminance(1.0)
                .values()
        }));
        // let xyzw = self.xyz_raw(&space.white, None).set_illuminance(1.0);
        let xyzw = self.xyz(&space.white, None).set_illuminance(1.0);
        let decomp = rgb2xyz.lu();
        // unwrap: only used with library color spaces
        let rgbw = decomp.solve(&xyzw.xyzn).unwrap();
        for (i, mut col) in rgb2xyz.column_iter_mut().enumerate() {
            col *= rgbw[i];
        }
        rgb2xyz
    }

    /// Calculates the RGB to XYZ matrix, for a particular color space.
    pub fn xyz2rgb(&self, rgbspace: RgbSpace) -> Matrix3<f64> {
        // unwrap: only used with library color spaces
        self.rgb2xyz(&rgbspace).try_inverse().unwrap()
    }
}

// JS-WASM Interface code
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl ObserverData {}

#[cfg(test)]
mod obs_test {

    use super::Observer;
    use crate::prelude::{StdIlluminant, CIE1931};
    use crate::rgbspace::RgbSpace;
    use approx::assert_ulps_eq;
    use strum::IntoEnumIterator as _;

    #[test]
    fn test_spectral_locus() {
        let [x, y] = CIE1931
            .spectral_locus_by_nm(CIE1931.spectral_locus_nm_min())
            .unwrap()
            .chromaticity();
        assert_ulps_eq!(x, 0.17411, epsilon = 1E-5);
        assert_ulps_eq!(y, 0.00496, epsilon = 1E-5);

        let [x, y] = CIE1931
            .spectral_locus_by_nm(CIE1931.spectral_locus_nm_max())
            .unwrap()
            .chromaticity();
        assert_ulps_eq!(x, 0.73469, epsilon = 1E-5);
        assert_ulps_eq!(y, 0.26531, epsilon = 1E-5);
    }

    #[test]
    fn test_spectral_locus_min_max() {
        let min = CIE1931.spectral_locus_index_min();
        println!("{min}");
        let max = CIE1931.spectral_locus_index_max();
        println!("{max}");
    }

    #[test]
    fn test_spectral_locus_to_rgb() {
        // FIXME: Once Observer implements `EnumIter`, make this test perform the conversion
        // for all observers.
        let observer = Observer::Std1931;
        let nm_min = observer.data().spectral_locus_nm_min();
        let nm_max = observer.data().spectral_locus_nm_max();

        for nm in nm_min..=nm_max {
            let xyz = observer.data().spectral_locus_by_nm(nm).unwrap();

            for rgbspace in RgbSpace::iter() {
                let rgb = xyz.rgb(Some(rgbspace));
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
        let got = CIE1931.rgb2xyz(&crate::rgbspace::RgbSpace::SRGB);
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
        let got = CIE1931.xyz2rgb(crate::rgbspace::RgbSpace::SRGB);
        approx::assert_ulps_eq!(want, got, epsilon = 3E-4);
    }

    #[test]
    fn test_xyz_std_illuminants() {
        use crate::xyz::XYZ;
        use nalgebra as na;
        let XYZ {
            xyz,
            xyzn,
            observer,
        } = CIE1931.xyz_cie_table(&StdIlluminant::D65, Some(100.0));
        approx::assert_ulps_eq!(xyzn, na::Vector3::new(95.04, 100.0, 108.86), epsilon = 1E-2);
        assert!(xyz.is_none());
    }

    #[test]
    fn test_planckian_locus() {
        // see https://www.waveformlighting.com/tech/calculate-cie-1931-xy-coordinates-from-cct
        // for test data (not clear what CMF domain they use)
        let xy = CIE1931.xyz_planckian_locus(3000.0).chromaticity();
        approx::assert_abs_diff_eq!(&xy.as_ref(), &[0.43693, 0.40407].as_ref(), epsilon = 2E-5);

        let xy = CIE1931.xyz_planckian_locus(6500.0).chromaticity();
        approx::assert_abs_diff_eq!(&xy.as_ref(), &[0.31352, 0.32363].as_ref(), epsilon = 6E-5);
    }

    #[test]
    fn test_xyz_from_illuminant_x_fn() {
        let xyz = CIE1931.xyz_from_std_illuminant_x_fn(&StdIlluminant::D65, |_v| 1.0);
        let d65xyz = CIE1931.xyz_d65().xyzn;
        approx::assert_ulps_eq!(
            xyz,
            crate::xyz::XYZ::from_vecs(d65xyz, Some(d65xyz), crate::observer::Observer::Std1931)
        );
    }

    #[test]
    fn test_xyz_of_sample_with_standard_illuminant() {
        use crate::prelude::{StdIlluminant::D65 as d65, XYZ};
        let xyz = CIE1931.xyz(&d65, Some(&crate::colorant::Colorant::white()));
        approx::assert_ulps_eq!(xyz, CIE1931.xyz_from_std_illuminant_x_fn(&d65, |_| 1.0));

        let xyz = CIE1931.xyz(&d65, Some(&crate::colorant::Colorant::black()));
        approx::assert_ulps_eq!(xyz, CIE1931.xyz_from_std_illuminant_x_fn(&d65, |_| 0.0));
    }
}
